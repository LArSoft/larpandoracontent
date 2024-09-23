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
#include <iterator>

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

    // Build up a map of every PFO to the slice it was reconstructed in.
    std::map<const ParticleFlowObject*, float> pfoToSliceIdMap;
    std::map<float, std::map<HitType, CaloHitList>> sliceIdToTotalHits;

    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        const PropertiesMap &properties(pPfo->GetPropertiesMap());
        const PropertiesMap::const_iterator iter(properties.find("SliceIndex"));

        if (iter != properties.end()) {
            pfoToSliceIdMap.insert({pPfo, iter->second});

            for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W}) {
                CaloHitList caloHits;
                LArPfoHelper::GetCaloHits(pPfo, view, caloHits);
                sliceIdToTotalHits[iter->second][view].insert(
                    sliceIdToTotalHits[iter->second][view].end(),
                    caloHits.begin(), caloHits.end()
                );
            }
        } else
            pfoToSliceIdMap.insert({pPfo, -1});

        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pRecoNeutrino = pPfo;
            break;
        }
    }

    const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};

    if (pRecoNeutrino && pTrueNeutrino)
    {
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRecoNeutrino)->GetPosition()};

        if (m_visualise)
        {
            const CartesianVector tu(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
            const CartesianVector tv(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
            const CartesianVector tw(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));

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

        if (m_writeFile && LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, m_detectorName))
        {
            const CartesianVector delta{recoVertex - trueVertex};
            const float dx{delta.GetX()}, dy{delta.GetY()}, dz{delta.GetZ()}, dr{delta.GetMagnitude()};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            const int success{1};

            // Was the true 3D vertex in any gap?
            const int isInGapUTrue(LArGeometryHelper::IsInGap3D(this->GetPandora(), trueVertex, TPC_VIEW_U));
            const int isInGapVTrue(LArGeometryHelper::IsInGap3D(this->GetPandora(), trueVertex, TPC_VIEW_V));
            const int isInGapWTrue(LArGeometryHelper::IsInGap3D(this->GetPandora(), trueVertex, TPC_VIEW_W));

            // Was the reco vertex, reconstructed in a gap?
            const int isInGapUReco(LArGeometryHelper::IsInGap3D(this->GetPandora(), recoVertex, TPC_VIEW_U));
            const int isInGapVReco(LArGeometryHelper::IsInGap3D(this->GetPandora(), recoVertex, TPC_VIEW_V));
            const int isInGapWReco(LArGeometryHelper::IsInGap3D(this->GetPandora(), recoVertex, TPC_VIEW_W));

            // Lets calculate some slice properties.
            const float sliceIndex(pfoToSliceIdMap.at(pRecoNeutrino));
            const auto mcHits = mcToHitsMap.at(pTrueNeutrino);
            CaloHitList uTrueHits, vTrueHits, wTrueHits;
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(uTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_U; });
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(vTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_V; });
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(wTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_W; });

            const auto uRecoHits = sliceIdToTotalHits.at(sliceIndex).at(TPC_VIEW_U);
            const auto vRecoHits = sliceIdToTotalHits.at(sliceIndex).at(TPC_VIEW_V);
            const auto wRecoHits = sliceIdToTotalHits.at(sliceIndex).at(TPC_VIEW_W);

            CaloHitList caloHits;
            LArPfoHelper::GetCaloHits(pRecoNeutrino, TPC_VIEW_U, caloHits);
            LArPfoHelper::GetCaloHits(pRecoNeutrino->GetDaughterPfoList(), TPC_VIEW_U, caloHits);
            const float nuCompU = LArMCParticleHelper::GetSharedHits(uTrueHits, caloHits).size() / (float) uTrueHits.size();
            const float nuPurityU = LArMCParticleHelper::GetSharedHits(uTrueHits, caloHits).size() / (float) uRecoHits.size();

            caloHits.clear();
            LArPfoHelper::GetCaloHits(pRecoNeutrino, TPC_VIEW_V, caloHits);
            LArPfoHelper::GetCaloHits(pRecoNeutrino->GetDaughterPfoList(), TPC_VIEW_V, caloHits);
            const float nuCompV = LArMCParticleHelper::GetSharedHits(vTrueHits, caloHits).size() / (float) vTrueHits.size();
            const float nuPurityV = LArMCParticleHelper::GetSharedHits(vTrueHits, caloHits).size() / (float) vRecoHits.size();

            caloHits.clear();
            LArPfoHelper::GetCaloHits(pRecoNeutrino, TPC_VIEW_W, caloHits);
            LArPfoHelper::GetCaloHits(pRecoNeutrino->GetDaughterPfoList(), TPC_VIEW_W, caloHits);
            const float nuCompW = LArMCParticleHelper::GetSharedHits(wTrueHits, caloHits).size() / (float) wTrueHits.size();
            const float nuPurityW = LArMCParticleHelper::GetSharedHits(wTrueHits, caloHits).size() / (float) wRecoHits.size();

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHits", mcHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoNuHitsU", uRecoHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoNuHitsV", vRecoHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoNuHitsW", wRecoHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapUTrue", isInGapUTrue));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapVTrue", isInGapVTrue));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapWTrue", isInGapWTrue));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapUReco", isInGapUReco));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapVReco", isInGapVReco));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapWReco", isInGapWReco));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSliceCompU", nuCompU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSlicePurU", nuPurityU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSliceCompV", nuCompV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSlicePurV", nuPurityV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSliceCompW", nuCompW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSlicePurW", nuPurityW));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
    else if (pTrueNeutrino)
    {
        if (m_writeFile && LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, m_detectorName))
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

StatusCode VertexMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DetectorName", m_detectorName));

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
