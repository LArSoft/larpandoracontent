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
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_filename;         // The filename of the ROOT output file
    std::string m_treename;         // The name of the ROOT tree
};

} // namespace lar_content

#endif // LAR_SLICE_MONITORING_ALGORITHM_H
