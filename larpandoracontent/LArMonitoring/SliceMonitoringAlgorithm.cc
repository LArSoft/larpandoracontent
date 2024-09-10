/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the slice monitoring algorithm.
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

void SliceMonitoringAlgorithm::RearrangeHits(const pandora::Algorithm *const pAlgorithm, SlicingAlgorithm::SliceList &inputSliceList,
    SlicingAlgorithm::SliceList &outputSliceList
){

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    CaloHitList fullCaloHitList{};
    for (const auto &slice : inputSliceList)
    {
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
            fullCaloHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
            fullCaloHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
            fullCaloHitList.push_back(pSliceCaloHit);
    }

    const MCParticle *pTrueNeutrino{nullptr};
    LArMCParticleHelper::CaloHitToMCMap caloHitToPrimaryMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;

    // Refactor: For slices find largest nu slice.
    //           For slices, populate the below variables.
    //             Can then add a "is best slice" tag for slice with most nu hits.
    //          Also add:
    //              IsBestSlice could be an int : 0 for yes, 1 for second best, 2
    //              for third... UniqueEventID : Lets me split up the slices a bit
    //              more, get an idea of how the event looks.
    //                  I.e. nuCompU where eventID == "ajK92kb" to see the
    //                  distribution of all the slices there.
    //              Double check the variables work for both "slice level" and "nu
    //              level"

    for (const CaloHit *pCaloHit : fullCaloHitList) {

        // Populate MC info.
        LArCaloHit *pLArCaloHit{
            const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};
        const auto mcWeights = pLArCaloHit->GetMCParticleWeightMap();

        const MCParticle *largestContributor{nullptr};
        float weight = -1;

        if (mcWeights.empty())
            continue;

        for (const auto &mcWeight : mcWeights) {
            const MCParticle *mc{mcWeight.first};
            const auto parent(LArMCParticleHelper::GetParentMCParticle(mc));

            if (mcWeight.second > weight)
                largestContributor = mc;

            if (LArMCParticleHelper::IsNeutrino(parent)) {
                pTrueNeutrino = parent;
                largestContributor = parent;
            }
        }

        if (largestContributor == nullptr)
            continue;

        caloHitToPrimaryMCMap[pCaloHit] = largestContributor;
        mcToTrueHitListMap[largestContributor].push_back(pCaloHit);
    }

    if (pTrueNeutrino) {
        const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
        const int success{1};

        // Lets calculate some slice properties.
        for (const auto &slice : inputSliceList)
        {
            CaloHitList sliceCaloHits;
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
                sliceCaloHits.push_back(pSliceCaloHit);
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
                sliceCaloHits.push_back(pSliceCaloHit);
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
                sliceCaloHits.push_back(pSliceCaloHit);

            CaloHitList mcHits;
            if (mcToTrueHitListMap.count(pTrueNeutrino) > 0)
                mcHits = mcToTrueHitListMap.at(pTrueNeutrino);

            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                                                   m_treename.c_str(), "trueNuHits",
                                                   (float)mcHits.size()));

            const std::map<HitType, std::string> allViews(
                {{TPC_VIEW_U, "_U"}, {TPC_VIEW_V, "_V"}, {TPC_VIEW_W, "_W"}});
            // std::map<std::string, std::vector<float>> scoreDistributions;

            for (const auto &viewNamePair : allViews) {
                HitType view(viewNamePair.first);
                auto viewName(viewNamePair.second);

                CaloHitList trueNuHitsForView;
                std::copy_if(mcHits.begin(), mcHits.end(),
                             std::back_inserter(trueNuHitsForView),
                             [&](const pandora::CaloHit *hit) {
                             return hit->GetHitType() == view;
                             });
                trueNuHitsForView.sort();

                CaloHitList allCaloHitsInView;
                std::copy_if(sliceCaloHits.begin(), sliceCaloHits.end(),
                             std::back_inserter(allCaloHitsInView),
                             [&](const pandora::CaloHit *hit) {
                             return hit->GetHitType() == view;
                             });
                allCaloHitsInView.sort();

                // Now we have all the hits + all the MC-based (i.e. neutrino hits), get
                // the overlap / missing / extra. First, the hits that match in the slice.
                CaloHitList sliceNuHits;
                std::set_intersection(trueNuHitsForView.begin(), trueNuHitsForView.end(),
                                      allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                      std::inserter(sliceNuHits, sliceNuHits.end()));

                // Now, the MC hits that are missing (i.e. the neutrino hits that are
                // missing).
                CaloHitList missingNuHits;
                std::set_difference(trueNuHitsForView.begin(), trueNuHitsForView.end(),
                                    allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                    std::inserter(missingNuHits, missingNuHits.end()));

                // Finally, the inverse, the hits that are in the slice that aren't
                // associated with the neutrino (cosmic contamination).
                CaloHitList sliceCRHits;
                std::set_difference(allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                    trueNuHitsForView.begin(), trueNuHitsForView.end(),
                                    std::inserter(sliceCRHits, sliceCRHits.end()));

                // Calculate slice completeness (how much of the neutrino is here), and
                // purity (how much CR contamination is there).
                const float containsNeutrinoHits = sliceNuHits.size() > 0;
                float nuComp(0.f);
                float nuPurity(1.f);

                if (containsNeutrinoHits) {
                    nuComp = sliceNuHits.size() / (float)trueNuHitsForView.size();
                    nuPurity = sliceNuHits.size() / (float)allCaloHitsInView.size();
                }

                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(),
                    "trueNuHitsInSlice" + viewName, (float)trueNuHitsForView.size()));
                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(), "sliceHits" + viewName,
                    (float)allCaloHitsInView.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "sliceNuHits" + viewName, (float)sliceNuHits.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "sliceCRHits" + viewName, (float)sliceCRHits.size()));
                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(), "missingNuHits" + viewName,
                    (float)missingNuHits.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "containsNeutrino" + viewName, containsNeutrinoHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                                                       m_treename.c_str(),
                                                       "nuSliceComp" + viewName, nuComp));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "nuSlicePur" + viewName, nuPurity));
            }
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    } else {
        // Lets calculate some slice properties.
        for (const auto &slice : inputSliceList)
        {
            const int success{0};
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    // Make the output slice list complete.
    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        outputSliceList.push_back(inputSliceList[sliceNumber]);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
