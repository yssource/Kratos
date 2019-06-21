//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez,
//


#include "project_sensitivity_to_skin_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "geometries/line_2d_2.h"



namespace Kratos
{
// Constructor for ProjectSensitivityToSkinProcess Process
ProjectSensitivityToSkinProcess::ProjectSensitivityToSkinProcess(ModelPart& rModelPart,
                    ModelPart& rSkinModelPart
                ):
    Process(),
    mrModelPart(rModelPart),
    mrSkinModelPart(rSkinModelPart)
{
}

void ProjectSensitivityToSkinProcess::Execute()
{
    KRATOS_TRY;

    BinBasedFastPointLocator<2> locator = BinBasedFastPointLocator<2>(mrModelPart);
    locator.UpdateSearchDatabase();
    const int max_results = 10000;
    typename BinBasedFastPointLocator<2>::ResultContainerType results(max_results);
    typename BinBasedFastPointLocator<2>::ResultIteratorType result_begin = results.begin();

    for(std::size_t i = 0; i < mrSkinModelPart.Nodes().size(); ++i) {
        auto it_skin_node=mrSkinModelPart.NodesBegin()+i;
        Vector NShapeFunc;
        Element::Pointer pElement;
        bool is_found = locator.FindPointOnMesh(it_skin_node->Coordinates(), NShapeFunc, pElement, result_begin, max_results);

        array_1d<double,3> normal = it_skin_node -> FastGetSolutionStepValue(NORMAL);
        double normal_norm = norm_2(normal);
        if (normal_norm < 1e-6){
            normal_norm = 1e-6;
        }
        array_1d<double,3> unit_normal = normal/normal_norm;

        if (is_found && pElement->Is(ACTIVE)){
            auto &r_geometry = pElement->GetGeometry();
            auto &r_sensitivity_vector = it_skin_node->FastGetSolutionStepValue(SHAPE_SENSITIVITY);
            r_sensitivity_vector.clear();

            for(std::size_t i_node = 0; i_node < 3; ++i_node) {
                double normal_sensitivity = r_geometry[i_node].FastGetSolutionStepValue(NORMAL_SENSITIVITY);
                for(std::size_t i_dim = 0; i_dim < 3; ++i_dim) {
                    r_sensitivity_vector[i_dim] += NShapeFunc[i_node]*normal_sensitivity*unit_normal[i_dim];
                }
            }
        }
    }


    KRATOS_CATCH("");
}
}// Namespace Kratos
