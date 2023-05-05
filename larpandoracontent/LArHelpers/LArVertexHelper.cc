/**
 *  @file   larpandoracontent/LArHelpers/LArVertexHelper.cc
 *
 *  @brief  Implementation of the vertex helper class.
 *
 *  $Log: $
 */

#include "Geometry/LArTPC.h"
#include "Managers/GeometryManager.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include <limits>

using namespace pandora;

namespace lar_content
{

LArVertexHelper::ClusterDirection LArVertexHelper::GetClusterDirectionInZ(const Pandora &pandora, const Vertex *const pVertex,
    const Cluster *const pCluster, const float tanAngle, const float apexShift)
{
    if ((VERTEX_3D != pVertex->GetVertexType()) || (tanAngle < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const CartesianVector theVertex2D(LArGeometryHelper::ProjectPosition(pandora, pVertex->GetPosition(), hitType));

    try
    {
        const LArPointingCluster pointingCluster(pCluster);
        const float length((pointingCluster.GetInnerVertex().GetPosition() - pointingCluster.GetOuterVertex().GetPosition()).GetMagnitude());
        const bool innerIsAtLowerZ(pointingCluster.GetInnerVertex().GetPosition().GetZ() < pointingCluster.GetOuterVertex().GetPosition().GetZ());

        float rLInner(std::numeric_limits<float>::max()), rTInner(std::numeric_limits<float>::max());
        float rLOuter(std::numeric_limits<float>::max()), rTOuter(std::numeric_limits<float>::max());
        LArPointingClusterHelper::GetImpactParameters(pointingCluster.GetInnerVertex(), theVertex2D, rLInner, rTInner);
        LArPointingClusterHelper::GetImpactParameters(pointingCluster.GetOuterVertex(), theVertex2D, rLOuter, rTOuter);

        const bool innerIsVertexAssociated(rLInner > (rTInner / tanAngle) - (length * apexShift));
        const bool outerIsVertexAssociated(rLOuter > (rTInner / tanAngle) - (length * apexShift));

        if (innerIsVertexAssociated == outerIsVertexAssociated)
            return DIRECTION_UNKNOWN;

        if ((innerIsVertexAssociated && innerIsAtLowerZ) || (outerIsVertexAssociated && !innerIsAtLowerZ))
            return DIRECTION_FORWARD_IN_Z;

        if ((innerIsVertexAssociated && !innerIsAtLowerZ) || (outerIsVertexAssociated && innerIsAtLowerZ))
            return DIRECTION_BACKWARD_IN_Z;
    }
    catch (StatusCodeException &)
    {
        return DIRECTION_UNKNOWN;
    }

    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsInFiducialVolume(const Pandora &pandora, const CartesianVector &vertex, const std::string &detector)
{
    const LArTPCMap &larTPCMap(pandora.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "VertexMonitoringAlgorithm::IsInFiducialVolume - LArTPC description not registered with Pandora as required " << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    float tpcMinX{std::numeric_limits<float>::max()}, tpcMaxX{-std::numeric_limits<float>::max()};
    float tpcMinY{std::numeric_limits<float>::max()}, tpcMaxY{-std::numeric_limits<float>::max()};
    float tpcMinZ{std::numeric_limits<float>::max()}, tpcMaxZ{-std::numeric_limits<float>::max()};

    for (const auto &[volumeId, pLArTPC] : larTPCMap)
    {
        (void)volumeId;
        const float centreX{pLArTPC->GetCenterX()}, halfWidthX{0.5f * pLArTPC->GetWidthX()};
        const float centreY{pLArTPC->GetCenterY()}, halfWidthY{0.5f * pLArTPC->GetWidthY()};
        const float centreZ{pLArTPC->GetCenterZ()}, halfWidthZ{0.5f * pLArTPC->GetWidthZ()};
        tpcMinX = std::min(tpcMinX, centreX - halfWidthX);
        tpcMaxX = std::max(tpcMaxX, centreX + halfWidthX);
        tpcMinY = std::min(tpcMinY, centreY - halfWidthY);
        tpcMaxY = std::max(tpcMaxY, centreY + halfWidthY);
        tpcMinZ = std::min(tpcMinZ, centreZ - halfWidthZ);
        tpcMaxZ = std::max(tpcMaxZ, centreZ + halfWidthZ);
    }

    const float x{vertex.GetX()};
    const float y{vertex.GetY()};
    const float z{vertex.GetZ()};

    if (detector == "dune_fd_hd")
    {
        return (tpcMinX + 50.f) < x && x < (tpcMaxX - 50.f) && (tpcMinY + 50.f) < y && y < (tpcMaxY - 50.f) && (tpcMinZ + 50.f) < z &&
               z < (tpcMaxZ - 150.f);
    }

    if (detector == "microboone")
    {
        return (tpcMinX + 10.f) < x && x < (tpcMaxX - 10.f) &&
               (tpcMinY + 10.f) < y && y < (tpcMaxY - 10.f) &&
               (tpcMinZ + 10.f) < z && z < (tpcMaxZ - 50.f);
    }

    if (detector == "microboone_dead_region")
    {
        return (tpcMinX + 10.f) < x && x < (tpcMaxX - 10.f) &&
               (tpcMinY + 10.f) < y && y < (tpcMaxY - 10.f) &&
               (tpcMinZ + 10.f) < z && z < (tpcMaxZ - 50.f) &&
               !((675 >= z) && (z <= 775)); // Additional restraint to ignore dead region.
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_content
