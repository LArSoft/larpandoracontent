/**
 *  @file   larpandoradlcontent/LArSliceTagging/DlEventSlicingTool.cc
 *
 *  @brief  Implementation of the deep learning slicing tool class.
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoradlcontent/LArSlicing/DlEventSlicingTool.h"
#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

typedef SlicingAlgorithm::HitTypeToNameMap HitTypeToNameMap;
typedef SlicingAlgorithm::SliceList SliceList;
typedef SlicingAlgorithm::Slice Slice;

DlEventSlicingTool::DlEventSlicingTool(): m_nClasses(3) {}

DlEventSlicingTool::~DlEventSlicingTool()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlEventSlicingTool: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlEventSlicingTool::RunSlicing(const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames, const HitTypeToNameMap &clusterListNames,
    SliceList &sliceList)
{
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};

    for (const auto &hitTypeListNamePair : caloHitListNames)
    {
        const auto listName = hitTypeListNamePair.second;
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_THROW_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, 
            PandoraContentApi::GetList(*pAlgorithm, listName, pCaloHitList)
        );

        if (pCaloHitList == nullptr || pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    Slice neutrinoSlice;

    for (const auto &hitTypeListNamePair : caloHitListNames)
    {
        HitType view{hitTypeListNamePair.first};
        const auto listName = hitTypeListNamePair.second;

        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_THROW_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
            PandoraContentApi::GetList(*pAlgorithm, listName, pCaloHitList)
        );

        if (pCaloHitList == nullptr || pCaloHitList->empty())
            continue;

        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return;

        LArDLHelper::TorchInput input;
        PixelVector pixelVector;
        CaloHitToPixelMap caloHitToPixelMap;
        this->MakeNetworkInputFromHits(*pCaloHitList, view, driftMin, driftMax, wireMin[view], wireMax[view], input, pixelVector, &caloHitToPixelMap);

        // Run the input through the trained model
        LArDLHelper::TorchInputVector inputs;
        inputs.push_back(input);
        LArDLHelper::TorchOutput output;

        if (isU)
            LArDLHelper::Forward(m_modelU, inputs, output);
        else if (isV)
            LArDLHelper::Forward(m_modelV, inputs, output);
        else
            LArDLHelper::Forward(m_modelW, inputs, output);

        // the argmax result is a 1 x height x width tensor where each element is a class id
        auto classesAccessor{output.accessor<float, 4>()};

        CaloHitList neutrinoSliceHits, otherHits;
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            auto found{caloHitToPixelMap.find(pCaloHit)};
            if (found == caloHitToPixelMap.end())
                continue;
            auto pixelMap = found->second;

            const int pixelZ(std::get<0>(pixelMap));
            const int pixelX(std::get<1>(pixelMap));

            // Apply softmax to loss to get actual probability
            float probNull = classesAccessor[0][0][pixelZ][pixelX];
            float probNeutrino = classesAccessor[0][1][pixelZ][pixelX];
            float probOther = classesAccessor[0][2][pixelZ][pixelX];

            if (probNeutrino > probOther && probNeutrino > probNull)
                neutrinoSliceHits.push_back(pCaloHit);
            else if (probOther > probNeutrino && probOther > probNull)
                otherHits.push_back(pCaloHit);

            float recipSum = 1.f / (probNeutrino + probOther);
            // Adjust probabilities to ignore null hits and update LArCaloHit
            probNeutrino *= recipSum;
            probOther *= recipSum;
        }

        if (m_visualise)
        {
            const std::string neutrinoListName("NeutrinoSlice_" + listName);
            const std::string otherListName("OtherSlice_" + listName);
            const std::string backgroundListName("BackgroundHits_" + listName);

            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &neutrinoSliceHits, neutrinoListName, BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &otherHits, otherListName, RED));
        }

        switch (view) {
            case TPC_VIEW_U:
                neutrinoSlice.m_caloHitListU = neutrinoSliceHits;
                break;
            case TPC_VIEW_V:
                neutrinoSlice.m_caloHitListV = neutrinoSliceHits;
                break;
            case TPC_VIEW_W:
                neutrinoSlice.m_caloHitListW = neutrinoSliceHits;
                break;
            default:
                break;
        }
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlEventSlicingTool::MakeNetworkInputFromHits(const CaloHitList &caloHits, const HitType view, const float xMin,
    const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector,
    CaloHitToPixelMap *caloHitToPixelMap) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float driftStep{0.5f};

    // Determine the bin edges
    std::vector<double> xBinEdges(m_width + 1);
    std::vector<double> zBinEdges(m_height + 1);
    xBinEdges[0] = xMin - 0.5f * driftStep;
    const double dx = ((xMax + 0.5f * driftStep) - xBinEdges[0]) / m_width;
    for (int i = 1; i < m_width + 1; ++i)
        xBinEdges[i] = xBinEdges[i - 1] + dx;
    zBinEdges[0] = zMin - 0.5f * pitch;
    const double dz = ((zMax + 0.5f * pitch) - zBinEdges[0]) / m_height;
    for (int i = 1; i < m_height + 1; ++i)
        zBinEdges[i] = zBinEdges[i - 1] + dz;

    LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    for (const CaloHit *pCaloHit : caloHits)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        const float adc{pCaloHit->GetMipEquivalentEnergy()};
        const int pixelX{static_cast<int>(std::floor((x - xBinEdges[0]) / dx))};
        const int pixelZ{static_cast<int>(std::floor((z - zBinEdges[0]) / dz))};
        accessor[0][0][pixelZ][pixelX] += adc;

        if (caloHitToPixelMap != nullptr)
            caloHitToPixelMap->insert({pCaloHit, {pixelZ, pixelX}});
    }
    for (int row = 0; row < m_height; ++row)
    {
        for (int col = 0; col < m_width; ++col)
        {
            const float value{accessor[0][0][row][col]};
            if (value > 0)
                pixelVector.emplace_back(std::make_pair(row, col));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlEventSlicingTool::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType view{caloHitList.front()->GetHitType()};
    const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
    if (!(isU || isV || isW))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    // Avoid unreasonable rescaling of very small hit regions, pixels are assumed to be 0.5cm in x and wire pitch in z
    // ATTN: Rescaling is to a size 1 pixel smaller than the intended image to ensure all hits fit within an imaged binned
    // to be one pixel wider than this
    const float xRange{xMax - xMin}, zRange{zMax - zMin};
    const float minXSpan{m_driftStep * (m_width - 1)};
    if (xRange < minXSpan)
    {
        const float padding{0.5f * (minXSpan - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    const float minZSpan{pitch * (m_height - 1)};
    if (zRange < minZSpan)
    {
        const float padding{0.5f * (minZSpan - zRange)};
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlEventSlicingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_height));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_width));

    std::string modelName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", modelName));
    modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
    LArDLHelper::LoadModel(modelName, m_modelU);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", modelName));
    modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
    LArDLHelper::LoadModel(modelName, m_modelV);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", modelName));
    modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
    LArDLHelper::LoadModel(modelName, m_modelW);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
