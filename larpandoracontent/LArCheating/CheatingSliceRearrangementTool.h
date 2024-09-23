/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceRearrangementTool.h
 *
 *  @brief  Header file for the cheating slice rearrangement tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SLICE_REARRANGEMENT_TOOL_H
#define LAR_CHEATING_SLICE_REARRANGEMENT_TOOL_H 1

#include "larpandoracontent/LArControlFlow/SlicingAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingSliceRearrangementTool class
 */
class CheatingSliceRearrangementTool : public SliceRearrangementBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingSliceRearrangementTool();

protected:
    /**
     *  @brief  Run the slicing tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputSliceList to receive the populated slice list
     *  @param  outputSliceList to receive the populated slice list
     */
    void RearrangeHits(const pandora::Algorithm *const pAlgorithm, SlicingAlgorithm::SliceList &inputSliceList, SlicingAlgorithm::SliceList &outputSliceList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:

    void AppendSliceHits(pandora::CaloHitList &hits, pandora::CaloHitList &hitsToFilter, pandora::HitType targetHitType);

    float m_targetSizeThreshold;   ///< The threshold to pick one slice as the "main" neutrino slice.
    float m_targetPurityThreshold; ///< Same but for "purity", i.e. percentage of nu hits to all hits.
    float m_moveSizeThreshold;     ///< The minimum target to consider moving hits into the "main" slice.
    float m_movePurityThreshold;   ///< Same but for "purity", i.e. percentage of nu hits to all hits.
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SLICE_REARRANGEMENT_TOOL_H
