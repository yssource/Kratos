//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

#if !defined(KRATOS_RANS_TEST_K_EPSILON_STABILIZED_CDR_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_TEST_K_EPSILON_STABILIZED_CDR_ADJOINT_UTILITIES_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"

#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"

#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_utilities.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_utilities.h"
#include "custom_utilities/test_stabilized_cdr_adjoint_interface_utilities.h"
#include "test_k_epsilon_utilities.h"

namespace Kratos
{
namespace Testing
{
template <typename TEvmElement, typename TEvmAdjointElement>
void RunScalarKEpsilon2D3NElementTest(
    const std::string PrimalElementName,
    const std::string AdjointElementName,
    const Variable<double>& rPerturbVariable,
    std::function<double(BoundedVector<double, 3>&,
                         const Vector&,
                         const Matrix&,
                         const typename TEvmAdjointElement::BaseType&,
                         const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValueAdjointSensitivities,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValue,
    const double Delta,
    const double RelativePrecision,
    const double AbsolutePrecision)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, PrimalElementName + "2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, AdjointElementName + "2D3N");

    ProcessInfo& r_adjoint_process_info = r_adjoint_model_part.GetProcessInfo();
    r_adjoint_process_info.SetValue(DELTA_TIME, r_adjoint_process_info[DELTA_TIME] * -1.0);

    std::vector<Process*> adjoint_processes_list;
    std::vector<Process*> primal_processes_list;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    adjoint_processes_list.push_back(&adjoint_nut_process);
    adjoint_processes_list.push_back(&nut_sensitivities_process);
    primal_processes_list.push_back(&primal_nut_process);

    auto perturbation = [rPerturbVariable](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(rPerturbVariable);
    };

    RunScalarSensitivityTest<TEvmElement, TEvmAdjointElement, 2, 3>(
        r_primal_model_part, r_adjoint_model_part, primal_processes_list,
        adjoint_processes_list, CalculateElementScalarValueAdjointSensitivities, perturbation,
        CalculateElementScalarValue, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        Delta, RelativePrecision, AbsolutePrecision);
}

template <typename TEvmElement, typename TEvmAdjointElement>
void RunVectorKEpsilon2D3NElementTest(
    const std::string PrimalElementName,
    const std::string AdjointElementName,
    const Variable<array_1d<double, 3>>& rPerturbVariable,
    std::function<double(BoundedMatrix<double, 3, 2>&,
                         const Vector&,
                         const Matrix&,
                         const typename TEvmAdjointElement::BaseType&,
                         const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValueAdjointSensitivities,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValue,
    const double Delta,
    const double RelativePrecision,
    const double AbsolutePrecision)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, PrimalElementName + "2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, AdjointElementName + "2D3N");

    ProcessInfo& r_adjoint_process_info = r_adjoint_model_part.GetProcessInfo();
    r_adjoint_process_info.SetValue(DELTA_TIME, r_adjoint_process_info[DELTA_TIME] * -1.0);

    std::vector<Process*> adjoint_processes_list;
    std::vector<Process*> primal_processes_list;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    adjoint_processes_list.push_back(&adjoint_nut_process);
    adjoint_processes_list.push_back(&nut_sensitivities_process);
    primal_processes_list.push_back(&primal_nut_process);

    auto perturbation = [rPerturbVariable](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_vector = rNode.FastGetSolutionStepValue(rPerturbVariable);
        return r_vector[Dim];
    };

    RunVectorSensitivityTest<TEvmElement, TEvmAdjointElement, 2, 3>(
        r_primal_model_part, r_adjoint_model_part, primal_processes_list,
        adjoint_processes_list, CalculateElementScalarValueAdjointSensitivities, perturbation,
        CalculateElementScalarValue, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        Delta, RelativePrecision, AbsolutePrecision);
}

template <typename TEvmElement, typename TEvmAdjointElement>
void RunShapeSensitivityKEpsilon2D3NElementTest(
    const std::string PrimalElementName,
    const std::string AdjointElementName,
    std::function<double(BoundedMatrix<double, 3, 2>&,
                         const Vector&,
                         const Matrix&,
                         const typename TEvmAdjointElement::BaseType&,
                         const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValueAdjointSensitivities,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValue,
    const double Delta,
    const double RelativePrecision,
    const double AbsolutePrecision)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, PrimalElementName + "2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, AdjointElementName + "2D3N");

    ProcessInfo& r_adjoint_process_info = r_adjoint_model_part.GetProcessInfo();
    r_adjoint_process_info.SetValue(DELTA_TIME, r_adjoint_process_info[DELTA_TIME] * -1.0);

    std::vector<Process*> adjoint_processes_list;
    std::vector<Process*> primal_processes_list;

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    adjoint_processes_list.push_back(&adjoint_nut_process);
    adjoint_processes_list.push_back(&nut_sensitivities_process);
    primal_processes_list.push_back(&primal_nut_process);

    auto perturbation = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    RunVectorSensitivityTest<TEvmElement, TEvmAdjointElement, 2, 3>(
        r_primal_model_part, r_adjoint_model_part, primal_processes_list,
        adjoint_processes_list, CalculateElementScalarValueAdjointSensitivities, perturbation,
        CalculateElementScalarValue, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        Delta, RelativePrecision, AbsolutePrecision);
}

} // namespace Testing
} // namespace Kratos

#endif