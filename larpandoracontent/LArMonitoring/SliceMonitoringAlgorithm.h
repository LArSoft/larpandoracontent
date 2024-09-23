/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h
 *
 *  @brief  Header file for the slice monitoring algorithm.
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
    void WriteOutHits(const std::map<unsigned int, pandora::CaloHitList> &inputSliceHits, const pandora::CaloHitList &threeDHits,
                      const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_trainingMode;              ///< Training mode
    std::string m_filename;           ///< The filename of the ROOT output file
    std::string m_treename;           ///< The name of the ROOT tree
    std::string m_trainingOutputFile; ///< Output name for training examples

    std::string m_trackPfoListName;   ///< Track PFOs to use for writing out training files
    std::string m_showerPfoListName;  ///< Shower PFOs to write out in the training files
};

} // namespace lar_content

#endif // LAR_SLICE_MONITORING_ALGORITHM_H
