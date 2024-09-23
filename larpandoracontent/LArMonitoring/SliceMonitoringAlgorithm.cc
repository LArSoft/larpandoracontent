/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the slice monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <iterator>

using namespace pandora;

namespace lar_content
{

SliceMonitoringAlgorithm::SliceMonitoringAlgorithm() : m_trainingMode(false), m_trackPfoListName(""), m_showerPfoListName("") {}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceMonitoringAlgorithm::~SliceMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceMonitoringAlgorithm::RearrangeHits(const pandora::Algorithm *const pAlgorithm, SlicingAlgorithm::SliceList &inputSliceList,
    SlicingAlgorithm::SliceList &outputSliceList
)
{

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    const CaloHitList *pCompleteCaloHitList{nullptr};
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetList(*pAlgorithm, "FullHitList", pCompleteCaloHitList));

    const PfoList *pTrackPfos(nullptr);
    if (PandoraContentApi::GetList(*pAlgorithm, m_trackPfoListName, pTrackPfos) != STATUS_CODE_SUCCESS)
        pTrackPfos = new PfoList();

    const PfoList *pShowerPfos(nullptr);
    if (PandoraContentApi::GetList(*pAlgorithm, m_showerPfoListName, pShowerPfos) != STATUS_CODE_SUCCESS)
        pShowerPfos = new PfoList();

    ClusterList *pClusterList(new ClusterList());
    for (const auto &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterListTemp(nullptr);
        PandoraContentApi::GetList(*pAlgorithm, clusterListName, pClusterListTemp);
        if (pClusterListTemp == nullptr)
            continue;
        pClusterList->insert(pClusterList->end(), pClusterListTemp->begin(), pClusterListTemp->end());
    }
    std::cout << "Cluster list has " << pClusterList->size() << " clusters." << std::endl;

    // Populate the complete calo hit list, based on every hit in every slice.
    CaloHitList fullSliceCaloHitList{};
    for (const auto &slice : inputSliceList)
    {
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
            fullSliceCaloHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
            fullSliceCaloHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
            fullSliceCaloHitList.push_back(pSliceCaloHit);
    }

    const MCParticle *pTrueNeutrino{nullptr};
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;

    // Also populate a list of 3D hits.
    // This is useful so we can assess the 3D completeness -> i.e. how much of a particle's hit were made into 3D hits.
    // Low 3D completeness == Hard to Slice
    CaloHitList threeDHits;
    for (const auto &pfoList : {pTrackPfos, pShowerPfos})
    {
        for (const auto pPfo : *pfoList )
        {
            CaloHitList pPfoHits;
            LArPfoHelper::GetCaloHits(pPfo, TPC_3D, pPfoHits);
            threeDHits.insert(threeDHits.end(), pPfoHits.begin(), pPfoHits.end());

            try
            {
                const auto pfoMC = LArMCParticleHelper::GetMainMCParticle(pPfo);
                for (const auto pCaloHit : pPfoHits)
                {
                    mcToTrueHitListMap[pfoMC].push_back(pCaloHit);
                }
            }
            catch (StatusCodeException &) { }

        }
    }

    // Since the slice hits, and the full, unfiltered hits may not match up, first
    // perform a quick match between these two sets of hits.
    typedef const std::tuple<HitType, float, float, float> CaloKey;
    std::map<CaloKey, const CaloHit*> caloHitListMatchMap;
    auto getHitKey = [](const CaloHit* pCaloHit) -> CaloKey {
        const auto pos = pCaloHit->GetPositionVector();
        return {pCaloHit->GetHitType(), pos.GetX(), pos.GetZ(), pCaloHit->GetHadronicEnergy()};
    };
    CaloHitList matchedCaloHitList;

    for (const CaloHit* pCaloHit : fullSliceCaloHitList)
        caloHitListMatchMap[getHitKey(pCaloHit)] = pCaloHit;

    for (const CaloHit *pCaloHit : *pCompleteCaloHitList) {

        // If there is a match for this hit...we want to point the slice-based hit,
        // to the complete calo hit list hit instead.
        // Otherwise, we can't really make connections between the two.
        if (caloHitListMatchMap.count(getHitKey(pCaloHit)) > 0)
            caloHitListMatchMap[getHitKey(pCaloHit)] = pCaloHit;

        // Populate MC info.
        LArCaloHit *pLArCaloHit{
            const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};
        const auto mcWeights = pLArCaloHit->GetMCParticleWeightMap();

        const MCParticle *largestContributor{nullptr};
        float weight = -1;

        if (mcWeights.empty())
            continue;

        for (const auto &mcWeight : mcWeights) {
            const MCParticle *mc{mcWeight.first};
            const auto parent(LArMCParticleHelper::GetParentMCParticle(mc));

            if (mcWeight.second > weight)
                largestContributor = mc;

            if (LArMCParticleHelper::IsNeutrino(parent)) {
                pTrueNeutrino = parent;
            }
        }

        if (largestContributor == nullptr)
            continue;

        mcToTrueHitListMap[largestContributor].push_back(pCaloHit);
    }

    if (pTrueNeutrino) {
        const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
        const int success{1};

        // WARN: This is valid at MicroBooNE, due to overlay information,
        //       but won't be valid if the cosmic rays have MC.
        CaloHitList mcHits;
        for (const auto &mcHitPair : mcToTrueHitListMap)
            mcHits.insert(mcHits.end(), mcHitPair.second.begin(), mcHitPair.second.end());
        // ATTN: Sort the list of MC hits, as the later std::set_XXX functions
        //       are undefined if run on unsorted sets.
        mcHits.sort();

        // First, assess the slices. We want to know which is the most-neutrino-filled, largest, etc.
        std::pair<unsigned int, unsigned int> bestSlice({0, 0});
        std::map<unsigned int, CaloHitList> matchedSliceHits;
        for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        {
            const auto slice(inputSliceList[sliceNumber]);
            CaloHitList sliceCaloHits({});

            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
                sliceCaloHits.push_back(caloHitListMatchMap[getHitKey(pSliceCaloHit)]);
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
                sliceCaloHits.push_back(caloHitListMatchMap[getHitKey(pSliceCaloHit)]);
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
                sliceCaloHits.push_back(caloHitListMatchMap[getHitKey(pSliceCaloHit)]);

            CaloHitList sliceNuHits;
            sliceCaloHits.sort();
            std::set_intersection(sliceCaloHits.begin(), sliceCaloHits.end(),
                                  mcHits.begin(), mcHits.end(),
                                  std::inserter(sliceNuHits, sliceNuHits.end()));

            if (sliceNuHits.size() > bestSlice.first)
                bestSlice = {sliceNuHits.size(), sliceNumber};

            matchedSliceHits[sliceNumber] = sliceCaloHits;
        }

        // INFO: Only want to write out training hits if we have a neutrino.
        if (m_trainingMode)
            WriteOutHits(matchedSliceHits, threeDHits, pClusterList, mcToTrueHitListMap);

        // Lets calculate some slice properties.
        // Perform this for every given input slice.
        for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        {
            const auto slice(inputSliceList[sliceNumber]);
            CaloHitList sliceCaloHits(matchedSliceHits[sliceNumber]);

            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "sliceNumber", (int) sliceNumber));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "mostNuHitSlice", (int) bestSlice.second));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "isBestSlice", (int) (sliceNumber == bestSlice.second)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                                                   m_treename.c_str(), "trueNuHits",
                                                   (float)mcHits.size()));

            const std::map<HitType, std::string> allViews(
                {{TPC_VIEW_U, "_U"}, {TPC_VIEW_V, "_V"}, {TPC_VIEW_W, "_W"}});

            for (const auto &viewNamePair : allViews) {
                HitType view(viewNamePair.first);
                auto viewName(viewNamePair.second);

                CaloHitList totalNuHitsInView;
                std::copy_if(mcHits.begin(), mcHits.end(),
                             std::back_inserter(totalNuHitsInView),
                             [&](const pandora::CaloHit *hit) {
                             return hit->GetHitType() == view;
                             });
                totalNuHitsInView.sort();

                CaloHitList allCaloHitsInView;
                std::copy_if(sliceCaloHits.begin(), sliceCaloHits.end(),
                             std::back_inserter(allCaloHitsInView),
                             [&](const pandora::CaloHit *hit) {
                             return hit->GetHitType() == view;
                             });
                allCaloHitsInView.sort();

                // Now we have all the hits + all the MC-based (i.e. neutrino hits), get
                // the overlap / missing / extra. First, the hits that match in the slice.
                CaloHitList sliceNuHits;
                std::set_intersection(totalNuHitsInView.begin(), totalNuHitsInView.end(),
                                      allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                      std::inserter(sliceNuHits, sliceNuHits.end()));

                // Now, the MC hits that are missing (i.e. the neutrino hits that are
                // missing).
                CaloHitList missingNuHits;
                std::set_difference(totalNuHitsInView.begin(), totalNuHitsInView.end(),
                                    allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                    std::inserter(missingNuHits, missingNuHits.end()));

                // Finally, the inverse, the hits that are in the slice that aren't
                // associated with the neutrino (cosmic contamination).
                CaloHitList sliceCRHits;
                std::set_difference(allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                    totalNuHitsInView.begin(), totalNuHitsInView.end(),
                                    std::inserter(sliceCRHits, sliceCRHits.end()));

                // Calculate slice completeness (how much of the neutrino is here), and
                // purity (how much CR contamination is there).
                const float containsNeutrinoHits = sliceNuHits.size() > 0;
                float nuComp(0.f);
                float nuPurity(1.f);

                if (containsNeutrinoHits) {
                    nuComp = sliceNuHits.size() / (float)totalNuHitsInView.size();
                    nuPurity = sliceNuHits.size() / (float)allCaloHitsInView.size();
                }

                std::cout << sliceNumber << ": " << nuComp << " / " << nuPurity <<
                        " (" << (sliceNumber == bestSlice.second) <<
                        " / " << (sliceNuHits.size()) <<
                        " / " << (sliceCRHits.size()) <<
                        " / " << (allCaloHitsInView.size()) <<
                        " / " << (missingNuHits.size()) <<
                        ")" << std::endl;

                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(),
                    "totalNuHitsInView" + viewName, (float)totalNuHitsInView.size()));
                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(), "sliceHits" + viewName,
                    (float)allCaloHitsInView.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "sliceNuHits" + viewName, (float)sliceNuHits.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "sliceCRHits" + viewName, (float)sliceCRHits.size()));
                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(), "missingNuHits" + viewName,
                    (float)missingNuHits.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "containsNeutrino" + viewName, containsNeutrinoHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                                                       m_treename.c_str(),
                                                       "nuSliceComp" + viewName, nuComp));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "nuSlicePur" + viewName, nuPurity));
            }

            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    } else {
        // Lets calculate some slice properties.
        for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        {
            const int success{0};
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "slice Number", (int) sliceNumber));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    // Make the output slice list complete.
    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        outputSliceList.push_back(inputSliceList[sliceNumber]);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceMonitoringAlgorithm::WriteOutHits(const std::map<unsigned int, CaloHitList> &inputSliceHits,
                                            const CaloHitList &threeDHits,
                                            const ClusterList *clusterList,
                                            const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap)
{

    // Build up a full list of all the slicing input hits, whilst retaining the
    // slice that Pandora currently chose.
    std::map<HitType, CaloHitList> perViewSliceHits({});
    std::map<const CaloHit*, int> hitToSliceNumber;

    for (const auto &sliceHitsPair : inputSliceHits )
    {
        const auto sliceNumber(sliceHitsPair.first);

        for (const CaloHit *const pSliceCaloHit : sliceHitsPair.second) {
            perViewSliceHits[pSliceCaloHit->GetHitType()].push_back(pSliceCaloHit);
            hitToSliceNumber.insert({pSliceCaloHit, sliceNumber});
        }
    }
    perViewSliceHits[TPC_3D] = threeDHits;

    std::map<const CaloHit*, int> hitToPdgCode;
    for (const auto &mcHitListPair : mcToTrueHitListMap)
        for (const auto &hit : mcHitListPair.second)
            hitToPdgCode.insert({hit, mcHitListPair.first->GetParticleId()});

    std::map<const std::tuple<float, float, float>, unsigned int> hitToClusterId;
    auto getHitKey = [](const CaloHit* pCaloHit) -> std::tuple<float, float, float> {
        const auto pos = pCaloHit->GetPositionVector();
        return {pos.GetX(), pos.GetZ(), pCaloHit->GetHadronicEnergy()};
    };
    unsigned int clusterId(1);
    for (auto it = clusterList->begin(); it != clusterList->end(); ++it, ++clusterId)
    {
        CaloHitList clusterHits;
        (*it)->GetOrderedCaloHitList().FillCaloHitList(clusterHits);

        for (const auto &hit : clusterHits)
            hitToClusterId.insert({getHitKey(hit), clusterId});
    }

    const std::map<HitType, std::string> allViews({
        {TPC_VIEW_U, "_U_View"}, {TPC_VIEW_V, "_V_View"},
        {TPC_VIEW_W, "_W_View"}, {TPC_3D, "_3D_Hits"}
    });

    for (const auto &viewNamePair : allViews)
    {
        HitType view(viewNamePair.first);
        auto viewName(viewNamePair.second);
        const auto sliceCaloHits(perViewSliceHits[view]);

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(inputSliceHits.size()));
        featureVector.emplace_back(static_cast<double>(sliceCaloHits.size()));

        for (const CaloHit *pCaloHit : sliceCaloHits)
        {
            const auto hitKey = getHitKey(pCaloHit);
            const float x{pCaloHit->GetPositionVector().GetX()};
            const float y{pCaloHit->GetPositionVector().GetY()};
            const float z{pCaloHit->GetPositionVector().GetZ()};
            const float adc{pCaloHit->GetMipEquivalentEnergy()};
            const float particlePdg(hitToPdgCode.count(pCaloHit) != 0 ? hitToPdgCode[pCaloHit] : 0);
            const float sliceNumber(hitToSliceNumber.count(pCaloHit) != 0 ? hitToSliceNumber[pCaloHit] : 0);
            const float clusterNum(hitToClusterId.count(hitKey) != 0 ? hitToClusterId[hitKey] : 0);

            featureVector.emplace_back(static_cast<double>(x));
            featureVector.emplace_back(static_cast<double>(y));
            featureVector.emplace_back(static_cast<double>(z));
            featureVector.emplace_back(static_cast<double>(adc));
            featureVector.emplace_back(static_cast<double>(particlePdg));
            featureVector.emplace_back(static_cast<double>(sliceNumber));
            featureVector.emplace_back(static_cast<double>(clusterNum));
        }

        const std::string trainingFileName{m_trainingOutputFile + viewName + ".csv"};
        LArMvaHelper::ProduceTrainingExample(trainingFileName, true, featureVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    if (m_trainingMode && (m_inputClusterListNames.empty()))
    {
        std::cout << "SliceMonitoring: Must provide cluster list names when producing training samples!" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));

    if (m_trainingMode)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
