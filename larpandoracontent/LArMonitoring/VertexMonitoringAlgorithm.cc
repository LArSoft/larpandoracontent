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
    std::cout << "Starting to assess vertices..." << std::endl;

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));

    // Get a mapping from each MC to its Calo hits...
    LArMCParticleHelper::CaloHitToMCMap caloHitToPrimaryMCMap;
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList2D, mcToPrimaryMCMap, caloHitToPrimaryMCMap, mcToHitsMap);

    CaloHitList mcHits;
    for (const auto &mcHitsPair : mcToHitsMap) {
        const auto hits = mcHitsPair.second;
        mcHits.insert(mcHits.end(), hits.begin(), hits.end());
    }

    // Then, find the primary MC particles.
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
                LArPfoHelper::GetCaloHits(pPfo->GetDaughterPfoList(), view, caloHits);
                sliceIdToTotalHits[iter->second][view].insert(
                    sliceIdToTotalHits[iter->second][view].end(),
                    caloHits.begin(), caloHits.end()
                );
            }
        }

        if (LArPfoHelper::IsNeutrino(pPfo) && pRecoNeutrino == nullptr)
        {
            pRecoNeutrino = pPfo;
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
            CaloHitList uTrueHits, vTrueHits, wTrueHits;
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(uTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_U; });
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(vTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_V; });
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(wTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_W; });

            const auto uSliceHits = sliceIdToTotalHits.at(sliceIndex).at(TPC_VIEW_U);
            const auto vSliceHits = sliceIdToTotalHits.at(sliceIndex).at(TPC_VIEW_V);
            const auto wSliceHits = sliceIdToTotalHits.at(sliceIndex).at(TPC_VIEW_W);

            CaloHitList caloHits;
            LArPfoHelper::GetCaloHits(pRecoNeutrino, TPC_VIEW_U, caloHits);
            LArPfoHelper::GetCaloHits(pRecoNeutrino->GetDaughterPfoList(), TPC_VIEW_U, caloHits);
            const int recoNuHitsU = caloHits.size();
            const float nuCompU = LArMCParticleHelper::GetSharedHits(uTrueHits, caloHits).size() / (float) uTrueHits.size();
            const float nuPurityU = LArMCParticleHelper::GetSharedHits(uTrueHits, caloHits).size() / (float) uSliceHits.size();
            const int isInSliceU = CheckIfSliceContainsVertex(uSliceHits, trueVertex, TPC_VIEW_U);

            caloHits.clear();
            LArPfoHelper::GetCaloHits(pRecoNeutrino, TPC_VIEW_V, caloHits);
            LArPfoHelper::GetCaloHits(pRecoNeutrino->GetDaughterPfoList(), TPC_VIEW_V, caloHits);
            const int recoNuHitsV = caloHits.size();
            const float nuCompV = LArMCParticleHelper::GetSharedHits(vTrueHits, caloHits).size() / (float) vTrueHits.size();
            const float nuPurityV = LArMCParticleHelper::GetSharedHits(vTrueHits, caloHits).size() / (float) vSliceHits.size();
            const int isInSliceV = CheckIfSliceContainsVertex(vSliceHits, trueVertex, TPC_VIEW_V);

            caloHits.clear();
            LArPfoHelper::GetCaloHits(pRecoNeutrino, TPC_VIEW_W, caloHits);
            LArPfoHelper::GetCaloHits(pRecoNeutrino->GetDaughterPfoList(), TPC_VIEW_W, caloHits);
            const int recoNuHitsW = caloHits.size();
            const float nuCompW = LArMCParticleHelper::GetSharedHits(wTrueHits, caloHits).size() / (float) wTrueHits.size();
            const float nuPurityW = LArMCParticleHelper::GetSharedHits(wTrueHits, caloHits).size() / (float) wSliceHits.size();
            const int isInSliceW = CheckIfSliceContainsVertex(wSliceHits, trueVertex, TPC_VIEW_W);

            // Top level Neutrino information
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHits", (int) mcHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsU", (int) uTrueHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsV", (int) vTrueHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsW", (int) wTrueHits.size()));

            // Then, vertex information.
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInSliceU", isInSliceU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInSliceV", isInSliceV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInSliceW", isInSliceW));

            // Finally, hit-level information. If the neutrino is badly reconstructed / missing, a worse vertex makes sense.
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoSliceHitsU", (int) uSliceHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoSliceHitsV", (int) vSliceHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoSliceHitsW", (int) wSliceHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoNuHitsU", recoNuHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoNuHitsV", recoNuHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoNuHitsW", recoNuHitsW));

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
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};

            // Was the true 3D vertex in any gap?
            const int isInGapUTrue(LArGeometryHelper::IsInGap3D(this->GetPandora(), trueVertex, TPC_VIEW_U));
            const int isInGapVTrue(LArGeometryHelper::IsInGap3D(this->GetPandora(), trueVertex, TPC_VIEW_V));
            const int isInGapWTrue(LArGeometryHelper::IsInGap3D(this->GetPandora(), trueVertex, TPC_VIEW_W));

            CaloHitList uTrueHits, vTrueHits, wTrueHits;
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(uTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_U; });
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(vTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_V; });
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(wTrueHits), [](const pandora::CaloHit* hit) { return hit->GetHitType() == TPC_VIEW_W; });

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHits", (int) mcHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsU", (int) uTrueHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsV", (int) vTrueHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsW", (int) wTrueHits.size()));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapUTrue", isInGapUTrue));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapVTrue", isInGapVTrue));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isInGapWTrue", isInGapWTrue));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheckIfSliceContainsVertex(const pandora::CaloHitList caloHits, const pandora::CartesianVector trueVertex, const HitType view) const
{
    const auto projectedVertex = LArGeometryHelper::ProjectPositionOntoView(trueVertex, view);

    float minX{std::numeric_limits<float>::max()};
    float maxX{std::numeric_limits<float>::min()};
    float minY{std::numeric_limits<float>::max()};
    float maxY{std::numeric_limits<float>::min()};
    float minZ{std::numeric_limits<float>::max()};
    float maxZ{std::numeric_limits<float>::min()};

    for (const auto &caloHit : caloHitList)
    {
        const auto pos = caloHit->GetPositionVector();
        minX = std::min(minX, pos.GetX());
        maxX = std::max(maxX, pos.GetX());
        minY = std::min(minY, pos.GetY());
        maxY = std::max(maxY, pos.GetY());
        minZ = std::min(minZ, pos.GetZ());
        maxZ = std::max(maxZ, pos.GetZ());
    }

    // TODO: Is just checking the x coordinate enough?
    if (minX == std::numeric_limits<float>::max() || maxX == std::numeric_limits<float>::min())
        return false;

    if (projectedVertex.GetX() < minX || projectedVertex.GetX() > maxX)
        return false;
    if (projectedVertex.GetY() < minY || projectedVertex.GetY() > maxY)
        return false;
    if (projectedVertex.GetZ() < minZ || projectedVertex.GetZ() > maxZ)
        return false;

    return true;
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
