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
#include "includes/cfd_variables.h"
#include "manufactured_solution_utility.h"
#include "manufactured_fluid_solutions_application_variables.h"
#include "processes/compute_nodal_gradient_process.h"
#include "processes/find_nodal_h_process.h"
#include "utilities/variable_utils.h"


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


void ManufacturedSolutionUtility::ComputeExactMaterialAcceleration()
{
    double time = mrModelPart.GetProcessInfo().GetValue(TIME);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        array_1d<double, 3>& coords = it_node->Coordinates();
        array_1d<double, 3> mat_acc = mrManufactured.TimeDerivative(coords, time) + mrManufactured.ConvectiveTerm(coords, time);
        it_node->SetValue(EXACT_MATERIAL_ACCELERATION, mat_acc);
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


void ManufacturedSolutionUtility::ComputeMaterialAccelerationError()
{
    ComputeRelativeError<Variable<array_1d<double, 3>>>(
        EXACT_MATERIAL_ACCELERATION,
        MATERIAL_ACCELERATION,
        MATERIAL_ACCELERATION_ERROR);
}


void ManufacturedSolutionUtility::RecoverMaterialAcceleration()
{
    VariableUtils().CopyVectorVar(ACCELERATION, MATERIAL_ACCELERATION, mrModelPart.Nodes());
    constexpr bool save_as_historical_variable = true;
    ComputeNodalGradientProcess<save_as_historical_variable>(mrModelPart, VELOCITY_X, VELOCITY_X_GRADIENT)();
    ComputeNodalGradientProcess<save_as_historical_variable>(mrModelPart, VELOCITY_Y, VELOCITY_Y_GRADIENT)();
    ComputeNodalGradientProcess<save_as_historical_variable>(mrModelPart, VELOCITY_Z, VELOCITY_Z_GRADIENT)();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        auto vel = it_node->FastGetSolutionStepValue(VELOCITY);
        auto acc = it_node->FastGetSolutionStepValue(ACCELERATION);
        auto grad_x = it_node->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
        auto grad_y = it_node->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
        auto grad_z = it_node->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
        it_node->FastGetSolutionStepValue(MATERIAL_ACCELERATION_X) += inner_prod(vel, grad_x);
        it_node->FastGetSolutionStepValue(MATERIAL_ACCELERATION_Y) += inner_prod(vel, grad_y);
        it_node->FastGetSolutionStepValue(MATERIAL_ACCELERATION_Z) += inner_prod(vel, grad_z);
    }
}


void ManufacturedSolutionUtility::ComputeNodalCFL()
{
    FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(mrModelPart).Execute();
    double time_step = mrModelPart.GetProcessInfo()[DELTA_TIME];
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin() + i;
        double h = it_node->GetValue(NODAL_H);
        double velocity = norm_2(it_node->FastGetSolutionStepValue(VELOCITY));
        it_node->SetValue(CFL_NUMBER, velocity * time_step / h);
    }
}


inline std::ostream& operator << (std::ostream& rOStream, const ManufacturedSolutionUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos.


