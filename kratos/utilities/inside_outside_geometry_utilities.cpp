//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes
#include <limits>
#include <algorithm>
#include <cmath>

// External includes
#include <delaunator-cpp/delaunator.hpp>

// Project includes
#include "utilities/inside_outside_geometry_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
namespace InsideOutsideGeometryUtilities
{
void GetLocalCoordinatesTriangle2D(
    const double LocationX,
    const double LocationY,
    const double Node1X,
    const double Node1Y,
    const double Node2X,
    const double Node2Y,
    const double Node3X,
    const double Node3Y,
    double& rLocalCoordinateX,
    double& rLocalCoordinateY
    )
{
    const double J_00 = Node2X - Node1X;
    const double J_01 = Node3X - Node1X;
    const double J_10 = Node2Y - Node1Y;
    const double J_11 = Node3Y - Node1Y;
    const double det_J = J_00 * J_11 - J_01 * J_10;

    // Compute eta and xi
    const double eta = (J_10*(Node1X - LocationX) + J_00 * (LocationY - Node1Y)) / det_J;
    const double xi = (J_11*(LocationX - Node1X) + J_01 * (Node1Y - LocationY)) / det_J;

    rLocalCoordinateX = xi;
    rLocalCoordinateY = eta;
}

/***********************************************************************************/
/***********************************************************************************/

bool IsInsideTriangle2D(
    const double LocationX,
    const double LocationY,
    const double Node1X,
    const double Node1Y,
    const double Node2X,
    const double Node2Y,
    const double Node3X,
    const double Node3Y,
    const double Tolerance
    )
{
    double local_coordinate_x = 0.0;
    double local_coordinate_y = 0.0;
    GetLocalCoordinatesTriangle2D(LocationX, LocationY, Node1X, Node1Y, Node2X, Node2Y, Node3X, Node3Y, local_coordinate_x, local_coordinate_y);

    if ((local_coordinate_x >= (0.0 - Tolerance)) && (local_coordinate_x <= (1.0 + Tolerance))) {
        if ((local_coordinate_y >= (0.0 - Tolerance)) && (local_coordinate_y <= (1.0 + Tolerance))) {
            if ((local_coordinate_x + local_coordinate_y) <= (1.0 + Tolerance)) {
                return true;
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool IsInside2D(
    const double LocationX,
    const double LocationY,
    const std::vector<std::vector<double>>& rOuterLoops,
    const std::vector<std::vector<double>>& rInnerLoops
    )
{
    for (std::size_t i = 0; i < rInnerLoops.size(); ++i) {
        delaunator::Delaunator delaunators(rInnerLoops[i]);

        for (std::size_t i = 0; i < delaunators.triangles.size(); i += 3) {
            const bool check_is_inside = IsInsideTriangle2D(
                LocationX,
                LocationY,
                rInnerLoops[i][2 * delaunators.triangles[i]],
                rInnerLoops[i][2 * delaunators.triangles[i] + 1],
                rInnerLoops[i][2 * delaunators.triangles[i + 1]],
                rInnerLoops[i][2 * delaunators.triangles[i + 1] + 1],
                rInnerLoops[i][2 * delaunators.triangles[i + 2]],
                rInnerLoops[i][2 * delaunators.triangles[i + 2] + 1]
            );
            if (check_is_inside)
                return false;
        }
    }
    for (std::size_t i = 0; i < rOuterLoops.size(); ++i) {
        delaunator::Delaunator delaunators(rOuterLoops[i]);

        for (std::size_t i = 0; i < delaunators.triangles.size(); i += 3) {
            const bool check_is_inside = IsInsideTriangle2D(
                LocationX,
                LocationY,
                rOuterLoops[i][2 * delaunators.triangles[i]],
                rOuterLoops[i][2 * delaunators.triangles[i] + 1],
                rOuterLoops[i][2 * delaunators.triangles[i + 1]],
                rOuterLoops[i][2 * delaunators.triangles[i + 1] + 1],
                rOuterLoops[i][2 * delaunators.triangles[i + 2]],
                rOuterLoops[i][2 * delaunators.triangles[i + 2] + 1]
            );
            if (check_is_inside)
                return true;
        }
    }

    return false;
}

} // namespace InsideOutsideGeometryUtilities
}  // namespace Kratos
