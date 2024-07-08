/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceRearrangementTool.cc
 *
 *  @brief  Implementation of the cheating slice rearrangement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingSliceRearrangementTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingSliceRearrangementTool::CheatingSliceRearrangementTool() :
    m_targetSizeThreshold{200}, m_targetPurityThreshold{0.75}, m_moveSizeThreshold{25}, m_movePurityThreshold{0.25}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CheatingSliceRearrangementTool::RearrangeHits(const pandora::Algorithm *const /*pAlgorithm*/,
                                                   SlicingAlgorithm::SliceList &inputSliceList, SlicingAlgorithm::SliceList &outputSliceList)
{

    // Before anything else, copy over the full slices, so we can quit out early if needed.
    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        outputSliceList.push_back(inputSliceList[sliceNumber]);

    int bestSlice(-1);
    std::map<int, std::pair<CaloHitList, float>> sliceMetrics;
    std::map<int, CaloHitList> sliceHits;

    // Get all the hits in each slice, and then find out if they are neutrino or cosmic-ray induced.
    // We can then calculate the number of neutrino hits and the percentage of nu hits out of all the hits in the slice.
    //
    // That info can be used to pick a slice as the "main" neutrino slice, and move neutrino hits from the other slices
    // into that slice.
    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
    {
        const auto slice(inputSliceList[sliceNumber]);
        CaloHitList localHitList{};
        CaloHitList localNuHitList{};

        // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
            localHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
            localHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
            localHitList.push_back(pSliceCaloHit);

        for (const CaloHit *const pCaloHit : localHitList)
        {
            const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

            // INFO: At MicroBooNE, if there is any MC at all...its a nu-hit.
            if (hitMCParticleWeightMap.empty())
                continue;

            localNuHitList.push_back(pCaloHit);
        }

        const unsigned int nuLikeHitCount(localNuHitList.size());
        const float slicePurity(nuLikeHitCount / (float) localHitList.size());
        sliceHits[sliceNumber] = localHitList;

        sliceMetrics.insert({sliceNumber, {localNuHitList, slicePurity}});

        std::cout << "Slice " << sliceNumber << ": " << localNuHitList.size() << " / " << localHitList.size() << " ( " << slicePurity << " % )" << std::endl;

        if (nuLikeHitCount < m_targetSizeThreshold || slicePurity < m_targetPurityThreshold)
            continue;

        if (
            bestSlice == -1 ||
            nuLikeHitCount > sliceMetrics[bestSlice].first.size() ||
            (nuLikeHitCount == sliceMetrics[bestSlice].first.size() && slicePurity > sliceMetrics[bestSlice].second)
        )
        {
            bestSlice = sliceNumber;
        }
    }

    if (bestSlice == -1)
        return;

    std::cout << "The best slice is slice " << bestSlice << std::endl;

    CaloHitList caloHitsToMove{};

    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
    {

        if (sliceNumber == (unsigned) bestSlice) continue;

        const auto sliceNuHits(sliceMetrics[sliceNumber].first);
        const auto slicePurity(sliceMetrics[sliceNumber].second);
        CaloHitList localHitList{};

        if (slicePurity < m_movePurityThreshold || sliceNuHits.size() < m_moveSizeThreshold)
            continue;

        caloHitsToMove.insert(caloHitsToMove.end(), sliceNuHits.begin(), sliceNuHits.end());
    }

    std::cout << "Moving " << caloHitsToMove.size() << " hits to the best slice!" << std::endl;

    if (caloHitsToMove.empty())
        return;

    for (unsigned int sliceNumber = 0; sliceNumber < outputSliceList.size(); ++sliceNumber)
    {
        // If this is the target slice, move hits into it.
        // Alternatively, do a quick check to remove hits from other slices.
        if (sliceNumber == (unsigned) bestSlice) {
            this->AppendSliceHits(outputSliceList[sliceNumber].m_caloHitListU, caloHitsToMove, TPC_VIEW_U);
            this->AppendSliceHits(outputSliceList[sliceNumber].m_caloHitListV, caloHitsToMove, TPC_VIEW_V);
            this->AppendSliceHits(outputSliceList[sliceNumber].m_caloHitListW, caloHitsToMove, TPC_VIEW_W);
        } else {

            std::map<pandora::HitType, CaloHitList> filteredHits;

            for (auto pCaloHit : sliceHits[sliceNumber]) {

                if (std::find(caloHitsToMove.begin(), caloHitsToMove.end(), pCaloHit) != caloHitsToMove.end())
                    continue;

                filteredHits[pCaloHit->GetHitType()].push_back(pCaloHit);
            }

            outputSliceList[sliceNumber].m_caloHitListU = filteredHits[TPC_VIEW_U];
            outputSliceList[sliceNumber].m_caloHitListV = filteredHits[TPC_VIEW_V];
            outputSliceList[sliceNumber].m_caloHitListW = filteredHits[TPC_VIEW_W];
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
void CheatingSliceRearrangementTool::AppendSliceHits(CaloHitList &hits, CaloHitList &hitsToFilter, pandora::HitType targetHitType)
{
    for (const auto &hit : hitsToFilter)
    {
        if (hit->GetHitType() == targetHitType)
            hits.push_back(hit);
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingSliceRearrangementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
