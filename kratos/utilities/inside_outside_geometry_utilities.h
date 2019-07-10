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

#if !defined(INSIDE_OUTSIDE_GEOMETRY_UTILITIES_H_INCLUDED)
#define INSIDE_OUTSIDE_GEOMETRY_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

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
namespace GeometryUtilities
{
    /**
     * @brief This method computes the local coordinates of  node inside a triangle (2D)
     * @note This is faster than create an auxiliar triangle, the code is repeated
     * @param LocationX The coordinate X of the point of study
     * @param LocationY The coordinate Y of the point of study
     * @param Node1X The coordinate X of the node 1
     * @param Node1Y The coordinate Y of the node 1
     * @param Node2X The coordinate X of the node 2
     * @param Node2Y The coordinate Y of the node 2
     * @param Node3X The coordinate X of the node 3
     * @param Node3Y The coordinate Y of the node 3
     * @param rLocalCoordinateX The local coordinate X of the point of study
     * @param rLocalCoordinateY The local coordinate Y of the point of study
     */
    void KRATOS_API(KRATOS_CORE) GetLocalCoordinatesTriangle2D(
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
        );
    /**
     * @brief This method checks if coordinates of a node are inside of a triangle (2D)
     * @note This is faster than create an auxiliar triangle, the code is repeated
     * @param LocationX The coordinate X of the point of study
     * @param LocationY The coordinate Y of the point of study
     * @param Node1X The coordinate X of the node 1
     * @param Node1Y The coordinate Y of the node 1
     * @param Node2X The coordinate X of the node 2
     * @param Node2Y The coordinate Y of the node 2
     * @param Node3X The coordinate X of the node 3
     * @param Node3Y The coordinate Y of the node 3
     * @param Tolerance The tolerance
     */
    bool KRATOS_API(KRATOS_CORE) IsInsideTriangle2D(
        const double LocationX,
        const double LocationY,
        const double Node1X,
        const double Node1Y,
        const double Node2X,
        const double Node2Y,
        const double Node3X,
        const double Node3Y,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        );

    /**
     * @brief This method creates a triangle mesh from two vectors of inner and outer loops and check is a node is inside
     * @param LocationX The coordinate X of the point of study
     * @param LocationY The coordinate Y of the point of study
     * @param rOuterLoops The outer loops definition
     * @param rOuterLoops The inner loops definition
     */
    bool KRATOS_API(KRATOS_CORE) IsInside2D(
        const double LocationX,
        const double LocationY,
        const std::vector<std::vector<double>>& rOuterLoops,
        const std::vector<std::vector<double>>& rInnerLoops
        );

}; // namespace GeometryUtilities
}  // namespace Kratos

#endif // INSIDE_OUTSIDE_GEOMETRY_UTILITIES_H_INCLUDED
