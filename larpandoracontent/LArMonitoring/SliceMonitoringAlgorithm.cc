/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <iterator>

using namespace pandora;

namespace lar_content
{

SliceMonitoringAlgorithm::SliceMonitoringAlgorithm() {}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceMonitoringAlgorithm::~SliceMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::Run()
{
    this->AssessSlice();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::AssessSlice() const
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    const CaloHitList *pCaloFullHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "FullHitList", pCaloFullHitList));

    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);

    const MCParticle *pTrueNeutrino{nullptr};

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

    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;

    // Just point every MC to the neutrino, CR don't have MC.
    for (const MCParticle *mc : *pMCParticleList)
        mcToPrimaryMCMap[mc] = pTrueNeutrino;

    LArMCParticleHelper::CaloHitToMCMap caloHitToPrimaryMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloFullHitList, mcToPrimaryMCMap, caloHitToPrimaryMCMap, mcToTrueHitListMap);

    if (pTrueNeutrino)
    {
        const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
        const int success{1};

        // Lets calculate some slice properties.
        CaloHitList mcHits;
        if (mcToTrueHitListMap.count(pTrueNeutrino) > 0)
            mcHits = mcToTrueHitListMap.at(pTrueNeutrino);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHits", (float) mcHits.size()));

        const std::map<HitType, std::string> allViews({{TPC_VIEW_U, "_U"}, {TPC_VIEW_V, "_V"}, {TPC_VIEW_W, "_W"}});
        std::map<std::string, std::vector<float>> scoreDistributions;

        for (const auto &viewNamePair : allViews)
        {
            HitType view(viewNamePair.first);
            auto viewName(viewNamePair.second);

            CaloHitList trueNuHitsForView;
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(trueNuHitsForView), [&](const pandora::CaloHit* hit) { return hit->GetHitType() == view; });
            trueNuHitsForView.sort();

            CaloHitList allCaloHitsInView;
            std::copy_if(pCaloHitList->begin(), pCaloHitList->end(), std::back_inserter(allCaloHitsInView), [&](const pandora::CaloHit* hit) { return hit->GetHitType() == view; });
            allCaloHitsInView.sort();

            // Now we have all the hits + all the MC-based (i.e. neutrino hits), get the overlap / missing / extra.
            // First, the hits that match in the slice.
            CaloHitList sliceNuHits;
            std::set_intersection(trueNuHitsForView.begin(), trueNuHitsForView.end(),
                                  allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                  std::inserter(sliceNuHits, sliceNuHits.end())
            );

            // Now, the MC hits that are missing (i.e. the neutrino hits that are missing).
            CaloHitList missingNuHits;
            std::set_difference(trueNuHitsForView.begin(), trueNuHitsForView.end(),
                                allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                std::inserter(missingNuHits, missingNuHits.end())
            );

            // Finally, the inverse, the hits that are in the slice that aren't associated with the neutrino (cosmic contamination).
            CaloHitList sliceCRHits;
            std::set_difference(allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                trueNuHitsForView.begin(), trueNuHitsForView.end(),
                                std::inserter(sliceCRHits, sliceCRHits.end())
            );

            // Calculate slice completeness (how much of the neutrino is here), and purity (how much CR contamination is there).
            const float containsNeutrinoHits = sliceNuHits.size() > 0;
            float nuComp(0.f);
            float nuPurity(1.f);

            if (containsNeutrinoHits) {
                nuComp = sliceNuHits.size() / (float) trueNuHitsForView.size();
                nuPurity = sliceNuHits.size() / (float) allCaloHitsInView.size();
            }

            std::vector<float> nuHitTagScores({});
            std::vector<float> crHitTagScores({});
            int tagCorrect(0);

            for (const CaloHit *pCaloHit : allCaloHitsInView) {

                LArCaloHit *pLArCaloHit{const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};

                if (caloHitToPrimaryMCMap.count(pCaloHit) == 0) {
                    crHitTagScores.push_back(pLArCaloHit->GetTrackProbability());
                    tagCorrect += pLArCaloHit->GetTrackProbability() > 0.5;
                } else {
                    nuHitTagScores.push_back(pLArCaloHit->GetShowerProbability());
                    tagCorrect += pLArCaloHit->GetShowerProbability() > 0.5;
                }
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsInSlice" + viewName, (float) trueNuHitsForView.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "sliceHits" + viewName, (float) allCaloHitsInView.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "sliceNuHits" + viewName, (float) sliceNuHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "sliceCRHits" + viewName, (float) sliceCRHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "missingNuHits" + viewName, (float) missingNuHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "containsNeutrino" + viewName, containsNeutrinoHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSliceComp" + viewName, nuComp));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSlicePur" + viewName, nuPurity));

            // Slice hit tagging efficiencies
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hitTagCorrectPct" + viewName, (float) tagCorrect / allCaloHitsInView.size()));

            scoreDistributions["trueNuHitScores" + viewName] = nuHitTagScores;
            scoreDistributions["trueCRHitScores" + viewName] = crHitTagScores;
        }

        for (auto &labelScorePair : scoreDistributions) {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), labelScorePair.first, &labelScorePair.second));
        }

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }
    else
    {
        const int success{0};
        const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
