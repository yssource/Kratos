//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/inside_outside_geometry_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(InsideOutsideGeometryUtilitiesIsInside2D, KratosCoreFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

    // First we create the nodes
    r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    std::vector<std::vector<double>> outer_loop(1), inner_loop(0);

    // The vector to fill
    std::vector<double> coordinates;

    for(auto& r_node : r_model_part.Nodes()) {
        coordinates.push_back(r_node.X());
        coordinates.push_back(r_node.Y());
    }

    outer_loop[0] = coordinates;

    double coor_x, coor_y;

    coor_x = 0.5;
    coor_y = 0.5;

    const bool test_1 = InsideOutsideGeometryUtilities::IsInside2D(coor_x, coor_y, outer_loop, inner_loop);
    KRATOS_CHECK(test_1);

    coor_x = 10.5;
    coor_y = 10.5;

    const bool test_2 = InsideOutsideGeometryUtilities::IsInside2D(coor_x, coor_y, outer_loop, inner_loop);
    KRATOS_CHECK(!test_2);
}

}   // namespace Testing
}  // namespace Kratos.
