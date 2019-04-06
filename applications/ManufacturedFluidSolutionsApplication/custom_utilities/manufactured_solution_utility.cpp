//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "manufactured_solution_utility.h"


namespace Kratos
{

ManufacturedSolutionUtility::ManufacturedSolutionUtility(
    ModelPart& rModelPart,
    ManufacturedSolution& rManufactured)
     : mrModelPart(rModelPart)
     , mrManufactured(rManufactured)
{}


void ManufacturedSolutionUtility::SetBodyForce()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        array_1d<double, 3>& coords = it_node->Coordinates();
        it_node->FastGetSolutionStepValue(BODY_FORCE) = mrManufactured.BodyForce(coords, time);
    }
}


void ManufacturedSolutionUtility::SetVelocity()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        array_1d<double, 3>& coords = it_node->Coordinates();
        it_node->FastGetSolutionStepValue(VELOCITY) = mrManufactured.Velocity(coords, time);
    }
}


void ManufacturedSolutionUtility::SetPressure()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        array_1d<double, 3>& coords = it_node->Coordinates();
        it_node->FastGetSolutionStepValue(PRESSURE) = mrManufactured.Pressure(coords, time);
    }
}


void ManufacturedSolutionUtility::ComputeExactVelocity()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        array_1d<double, 3>& coords = it_node->Coordinates();
        it_node->SetValue(EXACT_VELOCITY, mrManufactured.Velocity(coords, time));
    }
}


void ManufacturedSolutionUtility::ComputeExactPressure()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        array_1d<double, 3>& coords = it_node->Coordinates();
        it_node->SetValue(EXACT_PRESSURE, mrManufactured.Pressure(coords, time));
    }
}


void ManufacturedSolutionUtility::ComputeVelocityRelativeError()
{
    ComputeRelativeError<Variable<array_1d<double, 3>>>(EXACT_VELOCITY, VELOCITY, VELOCITY_RELATIVE_ERROR);
}


void ManufacturedSolutionUtility::ComputePressureRelativeError()
{
    ComputeRelativeError<Variable<double>>(EXACT_PRESSURE, PRESSURE, PRESSURE_RELATIVE_ERROR);
}

}  // namespace Kratos.


