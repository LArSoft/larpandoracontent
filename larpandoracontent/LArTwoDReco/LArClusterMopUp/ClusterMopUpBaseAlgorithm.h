/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/ClusterMopUpBaseAlgorithm.h
 *
 *  @brief  Header file for the cluster mop up algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_MOP_UP_BASE_ALGORITHM_H
#define LAR_CLUSTER_MOP_UP_BASE_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/MopUpBaseAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ClusterMopUpBaseAlgorithm class
 */

class ClusterMopUpBaseAlgorithm : public MopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterMopUpBaseAlgorithm();

protected:
    virtual pandora::StatusCode Run();

    /**
     *  @brief  Get the two dimensional clusters contained in the input pfo list, divided into three different views
     *
     *  @param  clusterListU to receive the list of clusters in the u view
     *  @param  clusterListV to receive the list of clusters in the v view
     *  @param  clusterListW to receive the list of clusters in the w view
     */
    virtual void GetPfoClusterLists(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Get the two dimensional clusters contained in the input remant/daughter cluster lists, divided into three different views
     *
     *  @param  clusterListU to receive the list of clusters in the u view
     *  @param  clusterListV to receive the list of clusters in the v view
     *  @param  clusterListW to receive the list of clusters in the w view
     */
    virtual void GetDaughterClusterLists(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Get the two dimensional clusters contained in an input cluster list, divided into three different views
     *
     *  @param  inputClusterList the input cluster list
     *  @param  availabilityFlag only clusters with matching availability will be considered
     *  @param  clusterListU to receive the list of clusters in the u view
     *  @param  clusterListV to receive the list of clusters in the v view
     *  @param  clusterListW to receive the list of clusters in the w view
     */
    virtual void GetClusterLists(const pandora::ClusterList &inputClusterList, const bool availabilityFlag,
        pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Cluster mop up for a single view. This function is responsible for instructing pandora to make cluster alterations
     *
     *  @param  pfoClusters the list of pfo clusters
     *  @param  remnantClusters the list of remnant clusters
     */
    virtual void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters) const = 0;

    typedef std::unordered_map<const pandora::Cluster *, float> AssociationDetails;
    typedef std::unordered_map<const pandora::Cluster *, AssociationDetails> ClusterAssociationMap;

    /**
     *  @brief  Make the cluster merges specified in the cluster association map, using list name information in the cluster list name map
     *
     *  @param  clusterAssociationMap the cluster association map
     */
    virtual void MakeClusterMerges(const ClusterAssociationMap &clusterAssociationMap) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames; ///< The list of pfo list names
    bool m_excludePfosContainingTracks;   ///< Whether to exclude any pfos containing clusters flagged as fixed tracks
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_MOP_UP_BASE_ALGORITHM_H
