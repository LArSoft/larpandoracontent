/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceRearrangementTool.h
 *
 *  @brief  Header file for the cheating slice rearrangement tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SLICE_REARRANGEMENT_TOOL_H
#define LAR_CHEATING_SLICE_REARRANGEMENT_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingSliceRearrangementTool class
 */
class CheatingSliceRearrangementTool : public SliceSelectionBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingSliceRearrangementTool();

    /**
     *  @brief  Select which slice(s) to use
     *
     *  @param  pAlgorithm the address of the master instance, used to access MCParticles when in training mode
     *  @param  inputSliceVector the initial slice vector
     *  @param  outputSliceVector the output slice vector
     */
    void SelectSlices(const pandora::Algorithm *const pAlgorithm, const SliceVector &inputSliceVector, SliceVector &outputSliceVector);

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

protected:
    float m_targetSizeThreshold;   ///< The threshold to pick one slice as the "main" neutrino slice.
    float m_targetPurityThreshold; ///< Same but for "purity", i.e. percentage of nu hits to all hits.
    float m_moveSizeThreshold;     ///< The minimum target to consider moving hits into the "main" slice.
    float m_movePurityThreshold;   ///< Same but for "purity", i.e. percentage of nu hits to all hits.
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SLICE_REARRANGEMENT_TOOL_H
