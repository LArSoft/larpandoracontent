/**
 *  @file   larpandoradlcontent/LArSlicing/DlEventSlicingTool.h
 *
 *  @brief  Header file for the deep learning slicing tool class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_EVENT_SLICING_TOOL_H
#define LAR_DL_EVENT_SLICING_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/SlicingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/EventSlicingTool.h"

#include "larpandoracontent/LArControlFlow/SlicingAlgorithm.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include <map>
#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DlEventSlicingTool class
 */
class DlEventSlicingTool : public EventSlicingTool
{
public:

    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;
    typedef std::map<const pandora::CaloHit *, std::tuple<int, int>> CaloHitToPixelMap;
    typedef std::pair<int, int> Pixel; // A Pixel is a row, column pair
    typedef std::vector<Pixel> PixelVector;

    /**
     *  @brief Default constructor
     */
    DlEventSlicingTool();

    virtual ~DlEventSlicingTool();

protected:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    void RunSlicing(const pandora::Algorithm *const pAlgorithm, const SlicingAlgorithm::HitTypeToNameMap &caloHitListNames,
        const SlicingAlgorithm::HitTypeToNameMap &clusterListNames, SlicingAlgorithm::SliceList &sliceList);

    void TagHits(const pandora::Algorithm *const pAlgorithm, const SlicingAlgorithm::HitTypeToNameMap &caloHitListNames,
        SlicingAlgorithm::Slice &neutrinoSlice);

    /**
     *  @brief  Populate a root true with vertex information.
     */
    void PopulateRootTree() const;

private:

    /*
     *  @brief  Create input for the network from a calo hit list
     *
     *  @param  caloHits The CaloHitList from which the input should be made
     *  @param  view The wire plane view
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  networkInput The TorchInput object to populate
     *  @param  pixelVector The output vector of populated pixels
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeNetworkInputFromHits(const pandora::CaloHitList &caloHits, const pandora::HitType view, const float xMin,
        const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector,
        CaloHitToPixelMap *caloHitToPixelMap = nullptr) const;

    /*
     *  @brief  Determine the physical bounds associated with a CaloHitList.
     *
     *  @param  caloHitList The CaloHitList for which bounds should be determined
     *  @param  xMin The output minimum x value
     *  @param  xMax The output maximum x value
     *  @param  zMin The output minimum z value
     *  @param  zMax The output maximum z value
     *
     *  @return The StatusCode resulting from the function
     */
    void GetHitRegion(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const;


    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    std::string m_inputVertexListName;        ///< Input vertex list name if 2nd pass
    std::string m_outputVertexListName;       ///< Output vertex list name
    LArDLHelper::TorchModel m_modelU;         ///< The model for the U view
    LArDLHelper::TorchModel m_modelV;         ///< The model for the V view
    LArDLHelper::TorchModel m_modelW;         ///< The model for the W view
    int m_event;                              ///< The current event number
    int m_nClasses;                           ///< The number of distance classes
    int m_height;                             ///< The height of the images
    int m_width;                              ///< The width of the images
    float m_driftStep;                        ///< The size of a pixel in the drift direction in cm (most relevant in pass 2)
    bool m_visualise;                         ///< Whether or not to visualise the candidate vertices
    bool m_writeTree;                         ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName;               ///< The ROOT tree name
    std::string m_rootFileName;               ///< The ROOT file name
    std::mt19937 m_rng;                       ///< The random number generator
    std::string m_detector;                   ///< Name of the current detector, for FV cuts
};

} // namespace lar_dl_content

#endif // LAR_DL_EVENT_SLICING_TOOL_H
