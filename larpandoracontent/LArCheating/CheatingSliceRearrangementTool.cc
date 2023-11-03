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

void CheatingSliceRearrangementTool::SelectSlices(const pandora::Algorithm *const /*pAlgorithm*/, const SliceVector &inputSliceVector, SliceVector &outputSliceVector)
{
    int bestSlice(-1);
    std::map<int, std::pair<CaloHitList, float>> sliceMetrics;

    // Get all the hits in each slice, and then find out if they are neutrino or cosmic-ray induced.
    // We can then calculate the number of neutrino hits and the percentage of nu hits out of all the hits in the slice.
    //
    // That info can be used to pick a slice as the "main" neutrino slice, and move neutrino hits from the other slices
    // into that slice.
    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceVector.size(); ++sliceNumber)
    {
        const auto sliceHits(inputSliceVector[sliceNumber]);
        CaloHitList localHitList{};
        CaloHitList localNuHitList{};

        // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
        for (const CaloHit *const pSliceCaloHit : sliceHits)
            localHitList.push_back(static_cast<const CaloHit *>(pSliceCaloHit->GetParentAddress()));

        for (const CaloHit *const pCaloHit : localHitList)
        {
            const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

            // INFO: At MicroBooNE, if there is any MC at all...its a nu-hit.
            if (hitMCParticleWeightMap.empty())
                continue;

            localNuHitList.push_back(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));
        }

        const unsigned int nuLikeHitCount(localNuHitList.size());
        const float slicePurity(localHitList.size() / (float) nuLikeHitCount);

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

    CaloHitList caloHitsToMove{};

    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceVector.size(); ++sliceNumber)
    {
        const auto sliceNuHits(sliceMetrics[sliceNumber].first);
        const auto slicePurity(sliceMetrics[sliceNumber].second);
        CaloHitList localHitList{};

        if (slicePurity < m_movePurityThreshold || sliceNuHits.size() < m_moveSizeThreshold)
            continue;

        caloHitsToMove.insert(caloHitsToMove.end(), sliceNuHits.begin(), sliceNuHits.end());
    }

    for (int sliceNumber = 0; sliceNumber < inputSliceVector.size(); ++sliceNumber)
    {
        auto sliceHits(inputSliceVector[sliceNumber]);

        // If this is the target slice, move hits into it.
        // Alternatively, do a quick check to remove hits from other slices.
        if (sliceNumber == bestSlice) {
            sliceHits.insert(sliceHits.end(), caloHitsToMove.begin(), caloHitsToMove.end());
        } else {
            CaloHitList filteredHits{};
            for (auto pCaloHit : sliceHits) {
                if (std::find(caloHitsToMove.begin(), caloHitsToMove.end(), pCaloHit) != caloHitsToMove.end())
                    continue;

                filteredHits.push_back(pCaloHit);
            }

            sliceHits = filteredHits;
        }

        outputSliceVector.push_back(sliceHits);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingSliceRearrangementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    // PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Threshold", m_threshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
