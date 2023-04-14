/**
 *  @file   larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

// Temporary includes for SCE correction
#include "TFile.h"

using namespace pandora;

namespace lar_content
{

VertexMonitoringAlgorithm::VertexMonitoringAlgorithm() :
    m_visualise{true},
    m_writeFile{false},
    m_detectorName{"dune_fd_hd"},
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexMonitoringAlgorithm::~VertexMonitoringAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::Run()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(
            this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));
    }

    this->AssessVertices();

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::AssessVertices() const
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    const ParticleFlowObject *pRecoNeutrino{nullptr};
    if (!primaries.empty())
    {
        for (const MCParticle *primary : primaries)
        {
            const MCParticleList &parents{primary->GetParentList()};
            if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            {
                pTrueNeutrino = parents.front();
                break;
            }
        }
    }

    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pRecoNeutrino = pPfo;
            break;
        }
    }

    // Load the SCE correction information.
    TFile inFile(m_pathToSCEFile.c_str(), "READ");

    if (!inFile.IsOpen())
    {
        std::cout << "VertexMonitoringAlgorithm: Could not find the space charge effect file!" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    TH3F* hDx = (TH3F*)inFile.Get("hDx");
    TH3F* hDy = (TH3F*)inFile.Get("hDy");
    TH3F* hDz = (TH3F*)inFile.Get("hDz");

    // Correct the truth information, first by transforming them to match the expected form in the SCE map.
    const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
    const CartesianVector transformedVertex(TransformVertex(trueVertex));

    // Then calculate any offset with the loaded maps.
    const CartesianVector positionOffset(GetPositionOffset(transformedVertex, hDx, hDy, hDz));

    // Finally, correct the true position with the calculated offsets.
    const CartesianVector trueVertexSCE(trueVertex.GetX() - positionOffset.GetX(),
                                        trueVertex.GetY() - positionOffset.GetY(),
                                        trueVertex.GetZ() - positionOffset.GetZ());

    if (pRecoNeutrino && pTrueNeutrino)
    {
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRecoNeutrino)->GetPosition()};

        if (m_visualise)
        {
            const CartesianVector tu(trueVertexSCE.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertexSCE.GetY(), trueVertexSCE.GetZ())));
            const CartesianVector tv(trueVertexSCE.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertexSCE.GetY(), trueVertexSCE.GetZ())));
            const CartesianVector tw(trueVertexSCE.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertexSCE.GetY(), trueVertexSCE.GetZ())));

            const CartesianVector ru(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(recoVertex.GetY(), recoVertex.GetZ())));
            const CartesianVector rv(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(recoVertex.GetY(), recoVertex.GetZ())));
            const CartesianVector rw(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(recoVertex.GetY(), recoVertex.GetZ())));

            const float du{(ru - tu).GetMagnitude()};
            const float dv{(rv - tv).GetMagnitude()};
            const float dw{(rw - tw).GetMagnitude()};

            std::cout << "delta(u, v, w): (" << du << ", " << dv << "," << dw << ")" << std::endl;

            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tu, "U true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tv, "V true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tw, "W true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &ru, "U reco vertex", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rv, "V reco vertex", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rw, "W reco vertex", RED, 2));
        }

        if (m_writeFile && IsInFiducialVolume(trueVertexSCE))
        {
            const CartesianVector delta{recoVertex - trueVertexSCE};
            const float dx{delta.GetX()}, dy{delta.GetY()}, dz{delta.GetZ()}, dr{delta.GetMagnitude()};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            const int success{1};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
    else if (pTrueNeutrino)
    {
        if (m_writeFile && IsInFiducialVolume(trueVertexSCE))
        {
            const int success{0};
            const float dx{-999.f}, dy{-999.f}, dz{-999.f}, dr{-999.f};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexMonitoringAlgorithm::TransformVertex(const CartesianVector& trueVertex) const
{
    // First, transform the position to match the SCE maps
    // Correct as of August 2019 -> TODO: Check if changed.
    const float newX(2.50f - (2.50f / 2.56f) * (trueVertex.GetX() / 100.0f));
    const float newY((2.50f / 2.33f) * ((trueVertex.GetY() / 100.0f) + 1.165f));
    const float newZ((10.0f / 10.37f) * (trueVertex.GetZ() / 100.0f));

    return CartesianVector(newX, newY, newZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexMonitoringAlgorithm::GetPositionOffset(const CartesianVector& transformedVertex, TH3F* dX, TH3F* dY, TH3F* dZ) const
{
    // Interpolate to get the positions offsets.
    const float xOffset(dX->Interpolate(transformedVertex.GetX(), transformedVertex.GetY(), transformedVertex.GetZ()));
    const float yOffset(dY->Interpolate(transformedVertex.GetX(), transformedVertex.GetY(), transformedVertex.GetZ()));
    const float zOffset(dZ->Interpolate(transformedVertex.GetX(), transformedVertex.GetY(), transformedVertex.GetZ()));

    return CartesianVector(xOffset, yOffset, zOffset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexMonitoringAlgorithm::IsInFiducialVolume(const CartesianVector &vertex) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "VertexMonitoringAlgorithm::IsInFiducialVolume - LArTPC description not registered with Pandora as required " << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    float tpcMinX{std::numeric_limits<float>::max()}, tpcMaxX{-std::numeric_limits<float>::max()};
    float tpcMinY{std::numeric_limits<float>::max()}, tpcMaxY{-std::numeric_limits<float>::max()};
    float tpcMinZ{std::numeric_limits<float>::max()}, tpcMaxZ{-std::numeric_limits<float>::max()};

    for (const auto &[volumeId, pLArTPC] : larTPCMap)
    {
        (void)volumeId;
        const float centreX{pLArTPC->GetCenterX()}, halfWidthX{0.5f * pLArTPC->GetWidthX()};
        const float centreY{pLArTPC->GetCenterY()}, halfWidthY{0.5f * pLArTPC->GetWidthY()};
        const float centreZ{pLArTPC->GetCenterZ()}, halfWidthZ{0.5f * pLArTPC->GetWidthZ()};
        tpcMinX = std::min(tpcMinX, centreX - halfWidthX);
        tpcMaxX = std::max(tpcMaxX, centreX + halfWidthX);
        tpcMinY = std::min(tpcMinY, centreY - halfWidthY);
        tpcMaxY = std::max(tpcMaxY, centreY + halfWidthY);
        tpcMinZ = std::min(tpcMinZ, centreZ - halfWidthZ);
        tpcMaxZ = std::max(tpcMaxZ, centreZ + halfWidthZ);
    }

    const float x{vertex.GetX()};
    const float y{vertex.GetY()};
    const float z{vertex.GetZ()};

    if (m_detectorName == "dune_fd_hd")
    {
        return (tpcMinX + 50.f) < x && x < (tpcMaxX - 50.f) && (tpcMinY + 50.f) < y && y < (tpcMaxY - 50.f) && (tpcMinZ + 50.f) < z &&
               z < (tpcMaxZ - 150.f);
    }

    if (m_detectorName == "microboone")
    {
        return (tpcMinX + 10.f) < x && x < (tpcMaxX - 10.f) &&
               (tpcMinY + 10.f) < y && y < (tpcMaxY - 10.f) &&
               (tpcMinZ + 10.f) < z && z < (tpcMaxZ - 50.f);
    }

    if (m_detectorName == "microboone_dead_region")
    {
        return (tpcMinX + 10.f) < x && x < (tpcMaxX - 10.f) &&
               (tpcMinY + 10.f) < y && y < (tpcMaxY - 10.f) &&
               (tpcMinZ + 10.f) < z && z < (tpcMaxZ - 50.f) &&
               !((675 >= z) && (z <= 775)); // Additional restraint to ignore dead region.
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DetectorName", m_detectorName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SCEFilePath", m_pathToSCEFile));

    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
