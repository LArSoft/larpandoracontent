/**
 *  @file   larpandoradlcontent/LArSliceTagging/DlSliceHitTagAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning slice hit tagging.
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoradlcontent/LArSliceTagging/DlSliceHitTagAlgorithm.h"
#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSliceHitTagAlgorithm::DlSliceHitTagAlgorithm(){}

DlSliceHitTagAlgorithm::~DlSliceHitTagAlgorithm()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlSliceHitTagAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSliceHitTagAlgorithm::Run()
{
    // 4 Output classes, NULL, Nu Track, Nu Shower, Cosmic.
    m_nClasses = 4;

    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlSliceHitTagAlgorithm::PrepareTrainingSample()
{
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    MCParticleList hierarchy;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));

    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

        const CaloHitList *caloU{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListU", caloU));
        const CaloHitList *caloV{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListV", caloV));
        const CaloHitList *caloW{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListW", caloW));

        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), caloU, "Calo U", GRAY));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), caloV, "Calo V", GRAY));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), caloW, "Calo W", GRAY));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    // Get boundaries for hits and make x dimension common
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!(isU || isV || isW))
            return STATUS_CODE_NOT_ALLOWED;

        std::map<const CaloHit*, int> caloHitToPDGMap;

        for (const MCParticle *mc : hierarchy)
            for (const CaloHit* hit : mcToHitsMap[mc])
                caloHitToPDGMap.insert({hit, mc->GetParticleId()});

        if (caloHitToPDGMap.empty())
            continue;

        const std::string trainingFilename{m_trainingOutputFile + "_" + listname + ".csv"};
        unsigned long nHits{0};
        const unsigned int nuance{LArMCParticleHelper::GetNuanceCode(hierarchy.front())};

        // Calo hits
        double xMin{driftMin}, xMax{driftMax}, zMin{wireMin[view]}, zMax{wireMax[view]};

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(nuance));
        // Retain the hit region
        featureVector.emplace_back(xMin);
        featureVector.emplace_back(xMax);
        featureVector.emplace_back(zMin);
        featureVector.emplace_back(zMax);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
            const float pdg(caloHitToPDGMap.count(pCaloHit) != 0 ? caloHitToPDGMap[pCaloHit] : 0);
            featureVector.emplace_back(static_cast<double>(x));
            featureVector.emplace_back(static_cast<double>(z));
            featureVector.emplace_back(static_cast<double>(adc));
            featureVector.emplace_back(static_cast<double>(pdg));
            ++nHits;
        }

        featureVector.insert(featureVector.begin() + 5, static_cast<double>(nHits));

        // Only write out the feature vector if there were enough hits in the region of interest
        if (nHits > 10)
            LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSliceHitTagAlgorithm::Infer()
{
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};

    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

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

        CaloHitList backgroundHits, neutrinoHits, cosmicRayHits;
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
            float probTrack = classesAccessor[0][1][pixelZ][pixelX];
            float probShower = classesAccessor[0][2][pixelZ][pixelX];
            float probCosmic = classesAccessor[0][3][pixelZ][pixelX];
            float probNeutrino = probTrack + probShower;

            if (probNeutrino > probCosmic && probNeutrino > probNull)
                neutrinoHits.push_back(pCaloHit);
            else if (probCosmic > probNeutrino && probCosmic > probNull)
                cosmicRayHits.push_back(pCaloHit);
            else
                backgroundHits.push_back(pCaloHit);

            float recipSum = 1.f / (probNeutrino + probCosmic);
            // Adjust probabilities to ignore null hits and update LArCaloHit
            probShower *= recipSum;
            probTrack *= recipSum;

            // TODO: Awful hack for proof of principle.
            //       Store Neutrino score in Shower, CR in track for later use.
            //       Need a smart way to persist the slice level metadata.
            LArCaloHit *pLArCaloHit{const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};
            pLArCaloHit->SetShowerProbability(probNeutrino);
            pLArCaloHit->SetTrackProbability(probCosmic);

            // caloHitToEvdHit[listName].at(pCaloHit)->addProperties({
            //     {"Neutrino-Like", probNeutrino},
            //     {"Cosmic-Like", probCosmic},
            //     {"isNeutrino", (float) probNeutrino > probCosmic},
            //     {"isCosmicRay", (float) probNeutrino < probCosmic}
            // });
        }

        if (m_visualise)
        {
            const std::string neutrinoListName("NeutrinoHits_" + listName);
            const std::string cosmicListName("CosmicHits_" + listName);
            const std::string backgroundListName("BackgroundHits_" + listName);

            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &neutrinoHits, neutrinoListName, BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &cosmicRayHits, cosmicListName, RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &backgroundHits, backgroundListName, BLACK));
        }
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSliceHitTagAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return DlVertexingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
