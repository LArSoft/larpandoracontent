/**
 *  @file   larpandoracontent/LArHelpers/LArHierarchyHelper.h
 *
 *  @brief  Header file for the lar hierarchy helper class.
 *
 *  $Log: $
 */
#ifndef LAR_HIERARCHY_HELPER_H
#define LAR_HIERARCHY_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

namespace lar_content
{

/**
 *  @brief  LArHierarchyHelper class
 */
class LArHierarchyHelper
{
public:
    /**
     *  @brief   MCHierarchy class
     */
    class MCHierarchy
    {
    public:
        /**
         *  @brief   ReconstructabilityCriteria class
         */
        class ReconstructabilityCriteria
        {
        public:
            /**
             *  @brief  Default constructor
             */
            ReconstructabilityCriteria();

            /**
             *  @brief  Copy constructor
             */
            ReconstructabilityCriteria(const ReconstructabilityCriteria &obj);

            /**
             *  @brief  Constructor
             *
             *  @param  minHits The total minimum number of hits for a particle to be considered reconstructable
             *  @param  minHitsForGoodView The number of hits within a view for a particle to be considered reconstructable
             *  @param  minGoodViews The minimum number of good views for a particle to be considered reconstructable
             *  @param removeNeutrons Whether to remove neutrons and downstream particles from consideration
             */
            ReconstructabilityCriteria(const unsigned int minHits, const unsigned int minHitsForGoodView, const unsigned int minGoodViews,
                const bool removeNeutrons);

            const unsigned int m_minHits;            ///< the minimum number of primary good Hits
            const unsigned int m_minHitsForGoodView; ///< the minimum number of Hits for a good view
            const unsigned int m_minGoodViews;       ///< the minimum number of primary good views
            const bool m_removeNeutrons;             ///< whether to remove neutrons and their downstream particles
        };

        class Node;
        typedef std::vector<const Node *> NodeVector;
        typedef std::list<const Node *> NodeList;

        /**
         *  @brief  Node class
         */
        class Node
        {
        public:
            /**
             *  @brief  Create a node with a primary MC particle
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  pMCParticle The primary MC particle with which this node should be created
             */
            Node(const MCHierarchy &hierarchy, const pandora::MCParticle *pMCParticle);

            /**
             *  @brief  Create a node from a list of MC particles
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  mcParticleList The MC particle list with which this node should be created
             *  @parasm caloHitList The CaloHit list with which this node should be created
             */
            Node(const MCHierarchy &hierarchy, const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList);

            /**
             *  @brief Destructor
             */
            virtual ~Node();

            /**
             *  @brief  Return whether or not this node should be considered reconstructable
             *
             *  @return true if reconstructable, false otherwise
             */
            bool IsReconstructable() const;

            /**
             *  @brief  Recursively fill the hierarchy based on the criteria established for this MCHierarchy
             *
             *  @param  pRoot The MC particle acting as the root of the current branch of the hierarchy
             *  @param foldToLeadingShower Whether or not we're folding back to the leading shower particle
             */
            void FillHierarchy(const pandora::MCParticle *pRoot, const bool foldToLeadingShower);

            /**
             *  @brief  Fill this node by folding all descendent particles to this node
             *
             *  @param  pRoot The MC particle acting as the root of the current branch of the hierarchy
             */
            void FillFlat(const pandora::MCParticle *pRoot);

            /**
             *  @brief  Return the vector of children for this node
             *
             *  @return The vector of children
             */
            const NodeVector &GetChildren() const;

            /**
             *  @brief  Retrieve the MC particles associated with this node
             *
             *  @return The MC particles associated with this node
             */
            const pandora::MCParticleList &GetMCParticles() const;

            /**
             *  @brief  Retrieve the CaloHits associated with this node
             *
             *  @return The list of CaloHits associated with this node
             */
            const pandora::CaloHitList &GetCaloHits() const;

            /**
             *  @brief  Retrieve the PDG code for the leading particle in this node
             *
             *  @return The PDG code for the leading particle in this node
             */
            int GetParticleId() const;

            /**
             *  @brief  Check if this is a particle induced by a neutrino interaction
             *
             *  @return Whether or not this is neutrino induced
             */
            bool IsNeutrinoInduced() const;

            /**
             *  @brief  Check if this is a test beam particle
             *
             *  @return Whether or not this is a test beam particle
             */
            bool IsTestBeamParticle() const;

            /**
             *  @brief  Check if this is a cosmic ray particle
             *
             *  @return Whether or not this is a cosmic ray
             */
            bool IsCosmicRay() const;

            /**
             *  @brief  Produce a string representation of the hierarchy
             *
             *  @return The string representation of the hierarchy
             */
            const std::string ToString(const std::string &prefix) const;

        private:
            const MCHierarchy &m_hierarchy;            ///< The parent MC hierarchy
            pandora::MCParticleList m_mcParticles;     ///< The list of MC particles of which this node is composed
            pandora::CaloHitList m_caloHits;           ///< The list of calo hits of which this node is composed
            NodeVector m_children;                     ///< The child nodes of this node
            const pandora::MCParticle *m_mainParticle; ///< The leading MC particle for this node
            int m_pdg;                                 ///< The PDG code of the leading MC particle for this node
        };

        /**
         *  @brief  Default constructor
         */
        MCHierarchy() = default;

        /**
         *  @brief  Construct a new MCHierarchy object using specified reconstructability criteria
         *
         *  @param  recoCriteria The reconstructability criteria to be applied
         */
        MCHierarchy(const ReconstructabilityCriteria &recoCriteria);

        /**
         *  @brief Destructor
         */
        virtual ~MCHierarchy();

        /**
         *  @brief  Creates an MC hierarchy representation. Without folding this will be a mirror image of the standard MCParticle
         *          relationships. However, with folding options selected the hierarchy structure will group together MC particles into
         *          nodes based on the folding requirements.
         *
         *          If only folding back to primaries, the hierarchy will be relatively flat, with a top-level neutrino or test beam
         *          particle, if appropriate, and then a set of leaf nodes, one for each primary particles also containing the MC particles
         *          (and corresponding hits) from daughter particles.
         *
         *          If only folding back to leading shower particles, the hierarchy will largely mirror the standard MCParticle hierarchy,
         *          but, when a shower particle is reached (for this purpose an electron or photon), this particle and all daughter
         *          particles will be represented by a single leaf node.
         *
         *          If folding back to both primary and leading shower particles the hierarchy will again be rather flat, but in this case,
         *          if a primary track-like particle (i.e. not an electron or photon) has a downstream shower particle then all downstream
         *          particles above the shower-like particle will be folded into the primary node, but a new, daughter leaf node will be
         *          created for the shower-like particle and all of its daughters, and a parent-child relationship will be formed between
         *          the primary node and shower node.
         *
         *  @param  mcParticleList The list of MC particles with which to fill the hierarchy
         *  @param  caloHitList The list of hits with which to fill the hierarchy
         *  @param  foldToPrimaries Whether or not to fold daughter particles back to their primary particle
         *  @param  foldToLeadingShowers Whether or not to fold daughter particles back to their leading shower particle
         */
        void FillHierarchy(const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList,
            const bool foldToPrimaries, const bool foldToLeadingShowers);

        /**
         *  @brief  Retrieve the root nodes in this hierarchy
         *
         *  @return The root nodes in this hierarchy
         */
        const NodeVector &GetRootNodes() const;

        /**
         *  @brief  Retrieve a flat vector of the ndoes in the hierarchy
         *
         *  @param  nodeVector The output vector for the nodes in the hierarchy in breadth first order
         */
        void GetFlattenedNodes(NodeVector &nodeVector) const;

        /**
         *  @brief  Produce a string representation of the hierarchy
         *
         *  @return The string representation of the hierarchy
         */
        const std::string ToString() const;

        /**
         *  @brief  Check if this is a neutrino hierarchy.
         *
         *  @return Whether or not this is a neutrino hierarchy.
         */
        bool IsNeutrinoHierarchy() const;

        /**
         *  @brief  Check if this is a test beam hierarchy.
         *
         *  @return Whether or not this is a test beam hierarchy.
         */
        bool IsTestBeamHierarchy() const;

    private:
        NodeVector m_rootNodes;                    ///< The leading nodes (e.g. primary particles, cosmic rays, ...)
        ReconstructabilityCriteria m_recoCriteria; ///< The criteria used to determine if the node is reconstructable
        const pandora::MCParticle *m_pNeutrino;    ///< The incident neutrino, if it exists
        std::map<const pandora::MCParticle *, pandora::CaloHitList> m_mcToHitsMap; ///< The map between MC particles and calo hits
    };

    /**
     *  @brief   RecoHierarchy class
     */
    class RecoHierarchy
    {
    public:
        class Node;
        typedef std::vector<const Node *> NodeVector;
        typedef std::list<const Node *> NodeList;

        /**
         *  @brief  Node class
         */
        class Node
        {
        public:
            /**
             *  @brief  Create a node with a primary PFO
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  pPfo The primary PFO with which this node should be created
             */
            Node(const RecoHierarchy &hierarchy, const pandora::ParticleFlowObject *pPfo);

            /**
             *  @brief  Create a node from a list of PFOs
             *
             *  @param  hierarchy The parent hierarchy of this node
             *  @param  pfoList The PFO list with which this node should be created
             *  @parasm caloHitList The CaloHit list with which this node should be created
             */
            Node(const RecoHierarchy &hierarchy, const pandora::PfoList &pfoList, const pandora::CaloHitList &caloHitList);

            /**
             *  @brief Destructor
             */
            virtual ~Node();

            /**
             *  @brief  Recursively fill the hierarchy based on the criteria established for this RecoHierarchy
             *
             *  @param  pRoot The PFO acting as the root of the current branch of the hierarchy
             *  @param foldToLeadingShower Whether or not we're folding back to the leading shower particle
             */
            void FillHierarchy(const pandora::ParticleFlowObject *pRoot, const bool foldToLeadingShower);

            /**
             *  @brief  Fill this node by folding all descendent particles to this node
             *
             *  @param  pRoot The PFO acting as the root of the current branch of the hierarchy
             */
            void FillFlat(const pandora::ParticleFlowObject *pRoot);

            /**
             *  @brief  Return the vector of children for this node
             *
             *  @return The vector of children
             */
            const NodeVector &GetChildren() const;

            /**
             *  @brief  Retrieve the PFOs associated with this node
             *
             *  @return The PFOs associated with this node
             */
            const pandora::PfoList &GetRecoParticles() const;

            /**
             *  @brief  Retrieve the CaloHits associated with this node
             *
             *  @return The list of CaloHits associated with this node
             */
            const pandora::CaloHitList &GetCaloHits() const;

            /**
             *  @brief  Retrieve the PDG code for the leading particle in this node
             *          Note, for reco objects the PDG codes represent tracks (muon PDG) and showers (electron PDG)
             *
             *  @return The PDG code for the leading particle in this node
             */
            int GetParticleId() const;

            /**
             *  @brief  Produce a string representation of the hierarchy
             *
             *  @return The string representation of the hierarchy
             */
            const std::string ToString(const std::string &prefix) const;

        private:
            const RecoHierarchy &m_hierarchy; ///< The parent reco hierarchy
            pandora::PfoList m_pfos;          ///< The list of PFOs of which this node is composed
            pandora::CaloHitList m_caloHits;  ///< The list of calo hits of which this node is composed
            NodeVector m_children;            ///< The child nodes of this node
            int m_pdg;                        ///< The particle ID (track = muon, shower = electron)
        };

        /**
         *  @brief  Default constructor
         */
        RecoHierarchy();

        /**
         *  @brief Destructor
         */
        virtual ~RecoHierarchy();

        /**
         *  @brief  Creates a reconstructed hierarchy representation. Without folding this will be a mirror image of the standard
         *          ParticleFlowObject (PFO) relationships. However, with folding options selected the hierarchy structure will group
         *          together PFOs into nodes based on the folding requirements.
         *
         *          If only folding back to primaries, the hierarchy will be relatively flat, with a top-level neutrino or test beam
         *          particle, if appropriate, and then a set of leaf nodes, one for each primary particles also containing the PFOs (and
         *          corresponding hits) from daughter particles.
         *
         *          If only folding back to leading shower particles, the hierarchy will largely mirror the standard PFO hierarchy, but,
         *          when a shower particle is reached (based on the track/shower characterisation), this particle and all daughter particles
         *          will be represented by a single leaf node.
         *
         *          If folding back to both primary and leading shower particles the hierarchy will again be rather flat, but in this case,
         *          if a primary track-like particle has a downstream shower particle then all downstream particles above the shower-like
         *          particle will be folded into the primary node, but a new, daughter leaf node will be created for the shower-like
         *          particle and all of its daughters, and a parent-child relationship will be formed between the primary node and shower
         *          node.
         *
         *  @param  pfoList The list of PFOs with which to fill the hierarchy
         *  @param  foldToPrimaries Whether or not to fold daughter particles back to their primary particle
         *  @param  foldToLeadingShowers Whether or not to fold daughter particles back to their leading shower particle
         */
        void FillHierarchy(const pandora::PfoList &pfoList, const bool foldToPrimaries, const bool foldToLeadingShowers);

        /**
         *  @brief  Retrieve the root nodes in this hierarchy
         *
         *  @return The root nodes in this hierarchy
         */
        const NodeVector &GetRootNodes() const;

        /**
         *  @brief  Retrieve a flat vector of the nodes in the hierarchy
         *
         *  @param  nodeVector The output vector for the nodes in the hierarchy in breadth first order
         */
        void GetFlattenedNodes(NodeVector &nodeVector) const;

        /**
         *  @brief  Produce a string representation of the hierarchy
         *
         *  @return The string representation of the hierarchy
         */
        const std::string ToString() const;

    private:
        NodeVector m_rootNodes;                         ///< The leading nodes (e.g. primary particles, cosmic rays, ...)
        const pandora::ParticleFlowObject *m_pNeutrino; ///< The incident neutrino, if it exists
    };

    /**
     *  @brief  MCMatches class
     */
    class MCMatches
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pMCParticle The MCParticle being matched
         */
        MCMatches(const MCHierarchy::Node *pMCParticle);

        /**
         *  @brief  Add a reconstructed node as a match for this MC node
         *
         *  @param  pReco The reconstructed node that matches this MC node
         *  @param  nSharedHits The number of hits shared betweeb reco and MC nodes
         */
        void AddRecoMatch(const RecoHierarchy::Node *pReco, const int nSharedHits);

        /**
         *  @brief  Retrieve the MC node
         *
         *  @return The MC node
         */
        const MCHierarchy::Node *GetMC() const;

        /**
         *  @brief  Retrieve the vector of matched reco nodes
         *
         *  @return The vector of matched reco nodes
         */
        const RecoHierarchy::NodeVector &GetRecoMatches() const;

        /**
         *  @brief  Retrieve the number of shared hits in the match
         *
         *  @param  pReco The reco node to consider
         *
         *  @return The number of shared hits
         */
        unsigned int GetSharedHits(const RecoHierarchy::Node *pReco) const;

        /**
         *  @brief  Retrieve the purity of the match
         *
         *  @param  pReco The reco node to consider
         *
         *  @return The purity of the match
         */
        float GetPurity(const RecoHierarchy::Node *pReco) const;

        /**
         *  @brief  Retrieve the completeness of the match
         *
         *  @param  pReco The reco node to consider
         *
         *  @return The completeness of the match
         */
        float GetCompleteness(const RecoHierarchy::Node *pReco) const;

    private:
        const MCHierarchy::Node *m_pMCParticle; ///< MC node associated with any matches
        RecoHierarchy::NodeVector m_recoNodes;  ///< Matched reco nodes
        pandora::IntVector m_sharedHits;        ///< Number of shared hits for each match
    };

    typedef std::vector<MCMatches> MCMatchesVector;

    /**
     *  @brief  MatcheInfo class
     */
    class MatchInfo
    {
    public:
        /**
         *  @brief  QualityCuts class
         */
        class QualityCuts
        {
        public:
            /**
             *  @brief Default constructor
             */
            QualityCuts();

            /**
             *  @brief Constructor
             *
             *  @param  minPurity The minimum purity for a cut to be considered good
             *  @param  minCompleteness The minimum completeness for a cut to be considered good
             */
            QualityCuts(const float minPurity, const float minCompleteness);

            const float m_minPurity;       ///< The minimum purity for a match to be considered good
            const float m_minCompleteness; ///< The minimum completeness for a match to be considered good
        };

        /**
         *  @brief  Default constructor
         */
        MatchInfo();

        /**
         *  @brief  Constructor
         *
         *  @param  qualityCuts The quality cuts to be applied to matched nodes
         */
        MatchInfo(const QualityCuts &qualityCuts);

        /**
         *  @brief  Match the nodes in the MC and reco hierarchies.
         *
         *  @param  mcHierarchy The MC hierarchy
         *  @param  recoHierarchy The reco hierarchy
         */
        void Match(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy);

        /**
         *  @brief  Retrieve the vector of good matches
         *
         *  @return The vector of good matches
         */
        const MCMatchesVector &GetGoodMatches() const;

        /**
         *  @brief  Retrieve the vector of matches that don't pass quality cuts
         *
         *  @return The vector of sub-threshold matches
         */
        const MCMatchesVector &GetSubThresholdMatches() const;

        /**
         *  @brief  Retrieve the vector of unmatched MC nodes
         *
         *  @return The vector of unmatched MC
         */
        const MCHierarchy::NodeVector &GetUnmatchedMC() const;

        /**
         *  @brief  Retrieve the vector of unmatched reco nodes
         */
        const RecoHierarchy::NodeVector &GetUnmatchedReco() const;

    private:
        MCMatchesVector m_goodMatches;             ///< The vector of good matches
        MCMatchesVector m_subThresholdMatches;     ///< The vector of matches MC nodes that don't pass quality cuts
        MCHierarchy::NodeVector m_unmatchedMC;     ///< The vector of unmatched MC nodes
        RecoHierarchy::NodeVector m_unmatchedReco; ///< The vector of unmatched reco nodes
        QualityCuts m_qualityCuts;                 ///< The quality cuts to be applied to matches
    };

    /**
     *  @brief  Fill an MC hierarchy based on the specified folding criteria (see MCHierarchy::FillHierarchy for details)
     *
     *  @param  mcParticleList The MCParticle list to use to fill this hierarchy
     *  @param  caloHitList The list of CaloHits to use to fill this hierarchy
     *  @param  foldToPrimaries Whether or not to fold to primary particles
     *  @param  foldToLeadingShowers Whether or not to fold to leading shower particles
     *  @param  hierarchy The output MC hierarchy
     */
    static void FillMCHierarchy(const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList,
        const bool foldToPrimaries, const bool foldToLeadingShowers, MCHierarchy &hierarchy);

    /**
     *  @brief  Fill a reconstructed hierarchy based on the specified folding criteria (see RecoHierarchy::FillHierarchy for details)
     *
     *  @param  pfoList The ParticleFlowObject list to use to fill this hierarchy
     *  @param  foldToPrimaries Whether or not to fold to primary particles
     *  @param  foldToLeadingShowers Whether or not to fold to leading shower particles
     *  @param  hierarchy The output reconstructed hierarchy
     */
    static void FillRecoHierarchy(const pandora::PfoList &pfoList, const bool foldToPrimaries, const bool foldToLeadingShowers, RecoHierarchy &hierarchy);

    /**
     *  @brief  Finds the matches between reconstructed and MC hierarchies.
     *
     *  @param  mcHierarchy The MC hiearchy
     *  @param  recoHierarchy The reconstructed hierarchy
     *  @param  matchInfo The output match information
     */
    static void MatchHierarchies(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy, MatchInfo &matchInfo);

private:
    typedef std::set<const pandora::MCParticle *> MCParticleSet;
    typedef std::set<const pandora::ParticleFlowObject *> PfoSet;

    /**
     *  @brief  Retrieves the primary MC particles from a list and returns the root (neutrino) for hierarchy, if it exists.
     *
     *  @param  mcParticleList The input list of MC particles
     *  @param  primaries The output set of primary MC particles
     *
     *  @return The root neutrino, if it exists, or nullptr
     */
    static const pandora::MCParticle *GetMCPrimaries(const pandora::MCParticleList &mcParticleList, MCParticleSet &primaries);

    /**
     *  @brief  Retrieves the primary PFOs from a list and returns the root (neutrino) for hierarchy, if it exists.
     *
     *  @param  pfoList The input list of PFOs
     *  @param  primaries The output set of primary PFOs
     *
     *  @return The root neutrino, if it exists, or nullptr
     */
    static const pandora::ParticleFlowObject *GetRecoPrimaries(const pandora::PfoList &pfoList, PfoSet &primaries);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy::NodeVector &LArHierarchyHelper::MCHierarchy::Node::GetChildren() const
{
    return m_children;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticleList &LArHierarchyHelper::MCHierarchy::Node::GetMCParticles() const
{
    return m_mcParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArHierarchyHelper::MCHierarchy::Node::GetCaloHits() const
{
    return m_caloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArHierarchyHelper::MCHierarchy::Node::GetParticleId() const
{
    return m_pdg;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArHierarchyHelper::MCHierarchy::Node::IsNeutrinoInduced() const
{
    return !(LArHierarchyHelper::MCHierarchy::Node::IsTestBeamParticle() || LArHierarchyHelper::MCHierarchy::Node::IsCosmicRay());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy::NodeVector &LArHierarchyHelper::MCHierarchy::GetRootNodes() const
{
    return m_rootNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArHierarchyHelper::MCHierarchy::IsNeutrinoHierarchy() const
{
    return m_pNeutrino != nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArHierarchyHelper::MCHierarchy::IsTestBeamHierarchy() const
{
    return m_pNeutrino == nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::RecoHierarchy::Node::GetChildren() const
{
    return m_children;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::RecoHierarchy::GetRootNodes() const
{
    return m_rootNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy::Node *LArHierarchyHelper::MCMatches::GetMC() const
{
    return m_pMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::MCMatches::GetRecoMatches() const
{
    return m_recoNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCMatchesVector &LArHierarchyHelper::MatchInfo::GetGoodMatches() const
{
    return m_goodMatches;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCMatchesVector &LArHierarchyHelper::MatchInfo::GetSubThresholdMatches() const
{
    return m_subThresholdMatches;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::MCHierarchy::NodeVector &LArHierarchyHelper::MatchInfo::GetUnmatchedMC() const
{
    return m_unmatchedMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArHierarchyHelper::RecoHierarchy::NodeVector &LArHierarchyHelper::MatchInfo::GetUnmatchedReco() const
{
    return m_unmatchedReco;
}

} // namespace lar_content

#endif // #ifndef LAR_HIERARCHY_HELPER_H
