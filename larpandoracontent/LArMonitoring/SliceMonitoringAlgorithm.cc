/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h"

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
    std::cout << "Assessing slice!" << std::endl;

    std::cout << "Getting MC..." << std::endl;
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    std::cout << "There is " << pMCParticleList->size() << " MC!" << std::endl;

    std::cout << "Getting hits..." << std::endl;
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    std::cout << "There is " << pCaloHitList->size() << " hits!" << std::endl;

    std::cout << "Populating maps..." << std::endl;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    std::cout << "Found " << primaries.size() << " primaries!" << std::endl;

    std::cout << "Finding neutrino..." << std::endl;
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
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToPrimaryMCMap, caloHitToPrimaryMCMap, mcToTrueHitListMap);
    std::cout << "There is " << mcToTrueHitListMap.size() << " MC entries..." << std::endl;
    std::cout << "There is " << caloHitToPrimaryMCMap.size() << " hit entries..." << std::endl;

    std::cout << "Writing stats out..." << std::endl;
    if (pTrueNeutrino)
    {
        const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
        const int success{1};

        // Lets calculate some slice properties.
        std::cout << "Getting the MC hits for the Nu..." << std::endl;
        const auto mcHits = mcToTrueHitListMap.at(pTrueNeutrino);

        std::cout << "Setting the basic properties..." << std::endl;
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHits", (float) mcHits.size()));

        const std::map<HitType, std::string> allViews({{TPC_VIEW_U, "_U"}, {TPC_VIEW_V, "_V"}, {TPC_VIEW_W, "_W"}});

        std::cout << "Starting loop..." << std::endl;
        for (const auto &viewNamePair : allViews)
        {
            HitType view(viewNamePair.first);
            auto viewName(viewNamePair.second);

            std::cout << viewName << std::endl;

            std::cout << "Getting MC Hits for view..." << std::endl;
            CaloHitList trueNuHitsForView;
            std::copy_if(mcHits.begin(), mcHits.end(), std::back_inserter(trueNuHitsForView), [&](const pandora::CaloHit* hit) { return hit->GetHitType() == view; });
            std::cout << "Got mc, there is " << trueNuHitsForView.size() << " MC hits..." << std::endl;
            trueNuHitsForView.sort();

            std::cout << "Getting Reco Hits for view..." << std::endl;
            CaloHitList allCaloHitsInView;
            std::copy_if(pCaloHitList->begin(), pCaloHitList->end(), std::back_inserter(allCaloHitsInView), [&](const pandora::CaloHit* hit) { return hit->GetHitType() == view; });
            std::cout << "Got hits, there is " << allCaloHitsInView.size() << " hits..." << std::endl;
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

            std::cout << "View" << viewName << ": " << nuComp << ", " << nuPurity << ", " << containsNeutrinoHits << std::endl;

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuHitsInSlice" + viewName, (float) trueNuHitsForView.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "sliceHits" + viewName, (float) allCaloHitsInView.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "sliceNuHits" + viewName, (float) sliceNuHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "sliceCRHits" + viewName, (float) sliceCRHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "missingNuHits" + viewName, (float) missingNuHits.size()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "containsNeutrino" + viewName, containsNeutrinoHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSliceComp" + viewName, nuComp));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuSlicePur" + viewName, nuPurity));

            // TODO: Hit tagging efficiencies
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
