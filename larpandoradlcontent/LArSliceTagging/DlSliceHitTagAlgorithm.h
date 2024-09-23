/**
 *  @file   larpandoradlcontent/LArSliceTagging/DlSliceHitTagAlgorithm.h
 *
 *  @brief  Header file for the deep learning slice hit tagging.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SLICE_HIT_TAG_ALGORITHM_H
#define LAR_DL_SLICE_HIT_TAG_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DlSliceHitTagAlgorithm class
 */
class DlSliceHitTagAlgorithm : public DlVertexingAlgorithm
{
public:

    /**
     *  @brief Default constructor
     */
    DlSliceHitTagAlgorithm();

    virtual ~DlSliceHitTagAlgorithm();

protected:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /**
     *  @brief  Populate a root true with vertex information.
     */
    void PopulateRootTree() const;
};

} // namespace lar_dl_content

#endif // LAR_DL_SLICE_HIT_TAG_ALGORITHM_H
