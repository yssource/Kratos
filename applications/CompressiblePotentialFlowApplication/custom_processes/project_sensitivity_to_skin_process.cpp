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

        if (is_found && pElement->Is(ACTIVE)){
            KRATOS_WATCH(NShapeFunc)
            KRATOS_WATCH(pElement->Id())
            KRATOS_WATCH(it_skin_node->Id())

            KRATOS_WATCH(it_skin_node->Coordinates())
        }
        // it_skin_node -> GetValue(SHAPE_SENSITIVITY) = it_skin_node->X();
    }

    for(std::size_t i = 0; i < mrSkinModelPart.Conditions().size(); ++i) {
        auto it_cond=mrSkinModelPart.ConditionsBegin()+i;
        auto r_geometry = it_cond -> GetGeometry();
        array_1d<double,3> aux_coords;
        // r_geometry.PrintInfo(std::cout);

        r_geometry.PointLocalCoordinates(aux_coords, r_geometry.Center());
        it_cond->SetValue(NORMAL, r_geometry.UnitNormal(aux_coords));
        // it_cond->SetValue(NORMAL, r_geometry.Normal());

    }


    KRATOS_CATCH("");
}
}// Namespace Kratos
