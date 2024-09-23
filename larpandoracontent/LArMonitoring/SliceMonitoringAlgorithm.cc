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

    const CaloHitList *pCaloSliceHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloSliceHitList));

    const CaloHitList *pCaloFullHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "FullHitList", pCaloFullHitList));

    const MCParticle *pTrueNeutrino{nullptr};
    LArMCParticleHelper::CaloHitToMCMap caloHitToPrimaryMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;

    // Since there is now two sets of calo hits...they will never match.
    // So build up a 2D map based on the position and find the matches
    // between the two lists manually.
    // First, populate the map, then use that in the next loop, looking
    // for exact matches.
    std::map<std::pair<float, float>, const CaloHit *> caloHitListMatchMap;
    CaloHitList pCaloHitListMatched;
    for (const CaloHit* pCaloHit : *pCaloSliceHitList)
    {
        const auto pos = pCaloHit->GetPositionVector();
        caloHitListMatchMap[{pos.GetX(), pos.GetZ()}] = pCaloHit;
    }

    for (const CaloHit* pCaloHit : *pCaloFullHitList)
    {
        // Now we can find the hits from the full list that match, and use them instead.
        const auto pos = pCaloHit->GetPositionVector();
        if (caloHitListMatchMap.count({pos.GetX(), pos.GetZ()}) > 0)
            if (caloHitListMatchMap[{pos.GetX(), pos.GetZ()}]->GetHadronicEnergy() == pCaloHit->GetHadronicEnergy())
                pCaloHitListMatched.push_back(pCaloHit);

        // Then populate MC info.
        LArCaloHit *pLArCaloHit{const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};
        const auto mcWeights = pLArCaloHit->GetMCParticleWeightMap();

        const MCParticle *largestContributor{nullptr};
        float weight;

        if (mcWeights.empty()) continue;

        for (const auto &mcWeight : mcWeights)
        {
            const MCParticle *mc{mcWeight.first};
            const auto parent(LArMCParticleHelper::GetParentMCParticle(mc));

            if (mcWeight.second > weight) largestContributor = mc;

            if (LArMCParticleHelper::IsNeutrino(parent))
            {
                pTrueNeutrino = parent;
                largestContributor = parent;
            }
        }

        if (largestContributor == nullptr) continue;

        caloHitToPrimaryMCMap[pCaloHit] = largestContributor;
        mcToTrueHitListMap[largestContributor].push_back(pCaloHit);
    }

    std::cout << "From " << pCaloFullHitList->size() << " hits, we ended up with " << caloHitToPrimaryMCMap.size() << " being mapped!" << std::endl;
    std::cout << "From " << pCaloSliceHitList->size() << " hits, we ended up with " << pCaloHitListMatched.size() << " being matched!" << std::endl;
    std::cout << "That is " << mcToTrueHitListMap.size() << " MCs matched!" << std::endl;

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
            std::copy_if(pCaloHitListMatched.begin(), pCaloHitListMatched.end(), std::back_inserter(allCaloHitsInView), [&](const pandora::CaloHit* hit) { return hit->GetHitType() == view; });
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

            if (containsNeutrinoHits)
            {
                nuComp = sliceNuHits.size() / (float) trueNuHitsForView.size();
                nuPurity = sliceNuHits.size() / (float) allCaloHitsInView.size();
            }

            std::vector<float> nuHitTagScores({});
            std::vector<float> crHitTagScores({});
            int tagCorrect(0);

            for (const CaloHit *pCaloHit : allCaloHitsInView)
            {
                const auto pos = pCaloHit->GetPositionVector();
                const CaloHit* sliceHit(caloHitListMatchMap[{pos.GetX(), pos.GetZ()}]);

                if (sliceHit->GetHadronicEnergy() != pCaloHit->GetHadronicEnergy())
                {
                    std::cout << "SliceVertexing: Slice to Full hit matching returned odd result!" << std::endl;
                    throw STATUS_CODE_INVALID_PARAMETER;
                }

                LArCaloHit *pLArCaloHit{const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(sliceHit))};

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
