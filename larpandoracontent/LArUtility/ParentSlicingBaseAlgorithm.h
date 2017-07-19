/**
 *  @file   larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h
 *
 *  @brief  Header file for the parent slicing base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_SLICING_BASE_ALGORITHM_H
#define LAR_PARENT_SLICING_BASE_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/ParentBaseAlgorithm.h"

namespace lar_content
{

class SlicingTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ParentSlicingBaseAlgorithm class
 */
class ParentSlicingBaseAlgorithm : public ParentBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ParentSlicingBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ParentSlicingBaseAlgorithm();

    /**
     *  @brief  Slice class
     */
    class Slice
    {
    public:
        pandora::CaloHitList    m_caloHitListU;                     ///< The u calo hit list
        pandora::CaloHitList    m_caloHitListV;                     ///< The v calo hit list
        pandora::CaloHitList    m_caloHitListW;                     ///< The w calo hit list
    };

    typedef std::vector<Slice> SliceList;

    /**
     *  @brief  Copy all the input hits in an event into a single slice
     *
     *  @param  sliceList the slice list to receive the single new slice
     */
    void CopyAllHitsToSingleSlice(SliceList &sliceList) const;

protected:
    /**
     *  @brief  Use first-pass 3D event reconstruction to slice events into separate, distinct interactions for processing
     *
     *  @param  sliceList the slice list to receive the slice list
     */
    void PerformSlicing(SliceList &sliceList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                        m_shouldPerformSlicing;             ///< Whether to slice events into separate, distinct interactions for processing
    SlicingTool                *m_pSlicingTool;                     ///< The address of the slicing tool
    std::string                 m_listDeletionAlgorithm;            ///< The name of the list deletion algorithm
    std::string                 m_listMovingAlgorithm;              ///< The name of the list moving algorithm
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SlicingTool class
 */
class SlicingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    virtual void Slice(const ParentSlicingBaseAlgorithm *const pAlgorithm, const ParentSlicingBaseAlgorithm::HitTypeToNameMap &caloHitListNames,
        const ParentSlicingBaseAlgorithm::HitTypeToNameMap &clusterListNames, ParentSlicingBaseAlgorithm::SliceList &sliceList) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_SLICING_BASE_ALGORITHM_H
