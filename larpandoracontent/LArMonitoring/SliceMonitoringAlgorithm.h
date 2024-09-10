/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_MONITORING_ALGORITHM_H
#define LAR_SLICE_MONITORING_ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/SlicingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  SliceMonitoringAlgorithm class
 */
class SliceMonitoringAlgorithm : public SliceRearrangementBaseTool
{
public:
    /**
   *  @brief  Default constructor
   */
    SliceMonitoringAlgorithm();

    virtual ~SliceMonitoringAlgorithm();

private:
    void RearrangeHits(const pandora::Algorithm *const pAlgorithm, SlicingAlgorithm::SliceList &inputSliceList, SlicingAlgorithm::SliceList &outputSliceList);
    void WriteOutHits(SlicingAlgorithm::SliceList &inputSliceList, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_trainingMode;              ///< Training mode
    std::string m_filename;           ///< The filename of the ROOT output file
    std::string m_treename;           ///< The name of the ROOT tree
    std::string m_trainingOutputFile; ///< Output name for training examples.
};

} // namespace lar_content

#endif // LAR_SLICE_MONITORING_ALGORITHM_H
