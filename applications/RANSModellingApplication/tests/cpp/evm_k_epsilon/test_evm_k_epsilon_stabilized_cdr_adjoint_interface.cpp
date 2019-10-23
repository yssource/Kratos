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
// System includes
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"

#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_element.h"

#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_utilities.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "test_k_epsilon_utilities.h"

#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"

namespace Kratos
{
namespace Testing
{
template <typename TEvmElement, typename TEvmAdjointElement, unsigned int TDim, unsigned int TNumNodes>
void RunScalarSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcessList,
    std::vector<Process*>& rAdjointProcessList,
    std::function<void(BoundedVector<double, TNumNodes>&,
                       const Vector&,
                       const Matrix&,
                       const typename TEvmAdjointElement::BaseType&,
                       const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                       const ProcessInfo&)> CalculateElementScalarValueAdjointSensitivities,
    std::function<double&(NodeType&)> PerturbVariable,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&)> CalculateElementScalarValue,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    const double Delta,
    const double RelativePrecision,
    const double AbsolutePrecision)
{
    const int number_of_primal_elements = rPrimalModelPart.NumberOfElements();
    const int number_of_adjoint_elements = rAdjointModelPart.NumberOfElements();

    KRATOS_ERROR_IF(number_of_primal_elements != number_of_adjoint_elements)
        << "Primal and Adjoint model part number of elements mismatch.\n";

    for (Process* p_process : rAdjointProcessList)
        p_process->Check();
    for (Process* p_process : rPrimalProcessList)
        p_process->Check();

    for (Process* p_process : rAdjointProcessList)
        p_process->ExecuteInitialize();
    for (Process* p_process : rAdjointProcessList)
        p_process->ExecuteInitializeSolutionStep();
    for (Process* p_process : rAdjointProcessList)
        p_process->Execute();

    UpdateVariablesInModelPart(rAdjointModelPart);

    auto execute_primal_process = [rPrimalProcessList]() {
        for (Process* p_process : rPrimalProcessList)
            p_process->ExecuteInitialize();

        for (Process* p_process : rPrimalProcessList)
            p_process->ExecuteInitializeSolutionStep();

        for (Process* p_process : rPrimalProcessList)
            p_process->Execute();
    };

    auto calculate_shape_functions =
        [](typename Element::GeometryType::ShapeFunctionsGradientsType& rShapeFunctionDerivatives,
           Element& rElement) {
            auto& r_geometry = rElement.GetGeometry();
            const auto& r_integration_method = rElement.GetIntegrationMethod();
            r_geometry.ShapeFunctionsIntegrationPointsGradients(
                rShapeFunctionDerivatives, r_integration_method);
            return r_geometry.ShapeFunctionsValues(r_integration_method);
        };

    const ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();
    const ProcessInfo& r_adjoint_process_info = rAdjointModelPart.GetProcessInfo();

    typename Element::GeometryType::ShapeFunctionsGradientsType adjoint_shape_function_gradients;
    typename Element::GeometryType::ShapeFunctionsGradientsType primal_shape_function_gradients;

    for (int i = 0; i < number_of_primal_elements; ++i)
    {
        ModelPart::ElementType& r_primal_element =
            *(rPrimalModelPart.ElementsBegin() + i);
        ModelPart::ElementType& r_adjoint_element =
            *(rAdjointModelPart.ElementsBegin() + i);

        const auto& r_primal_integration_method = r_primal_element.GetIntegrationMethod();
        const auto& r_adjoint_integration_method =
            r_adjoint_element.GetIntegrationMethod();

        KRATOS_ERROR_IF(r_primal_integration_method != r_adjoint_integration_method)
            << "Primal and Adjoint integration method mismatch.\n";

        const Matrix& adjoint_shape_functions = calculate_shape_functions(
            adjoint_shape_function_gradients, r_adjoint_element);

        auto& r_primal_geometry = r_primal_element.GetGeometry();

        auto& r_rans_adjoint_element =
            dynamic_cast<const typename TEvmAdjointElement::BaseType&>(r_adjoint_element);
        auto& r_rans_primal_element =
            dynamic_cast<const typename TEvmElement::BaseType&>(r_primal_element);

        const int number_of_gauss_points = adjoint_shape_functions.size1();
        for (int g = 0; g < number_of_gauss_points; ++g)
        {
            const Vector& r_adjoint_gauss_shape_functions =
                row(adjoint_shape_functions, g);
            const Matrix& r_adjoint_gauss_shape_derivatives =
                adjoint_shape_function_gradients[g];

            typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType adjoint_data;
            r_rans_adjoint_element.CalculateElementData(
                adjoint_data, r_adjoint_gauss_shape_functions,
                r_adjoint_gauss_shape_derivatives, r_adjoint_process_info);

            BoundedVector<double, TNumNodes> adjoint_scalar_sensitivities;
            CalculateElementScalarValueAdjointSensitivities(
                adjoint_scalar_sensitivities, r_adjoint_gauss_shape_functions,
                r_adjoint_gauss_shape_derivatives, r_rans_adjoint_element,
                adjoint_data, r_adjoint_process_info);

            // Calculating reference scalar value
            execute_primal_process();
            UpdateVariablesInModelPart(rPrimalModelPart);
            const Matrix& primal_shape_functions = calculate_shape_functions(
                primal_shape_function_gradients, r_primal_element);
            const Vector& primal_gauss_shape_functions = row(primal_shape_functions, g);
            const Matrix& primal_gauss_shape_derivatives =
                primal_shape_function_gradients[g];
            typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType data;
            r_rans_primal_element.CalculateConvectionDiffusionReactionData(
                data, primal_gauss_shape_functions,
                primal_gauss_shape_derivatives, r_primal_process_info);
            const double scalar_value_reference = CalculateElementScalarValue(
                primal_gauss_shape_functions, primal_gauss_shape_derivatives,
                r_rans_primal_element, data, r_primal_process_info);

            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            {
                NodeType& r_node = r_primal_geometry[i_node];
                PerturbVariable(r_node) += Delta;

                execute_primal_process();
                UpdateVariablesInModelPart(rPrimalModelPart);
                const Matrix& primal_shape_functions = calculate_shape_functions(
                    primal_shape_function_gradients, r_primal_element);
                const Vector& primal_gauss_shape_functions =
                    row(primal_shape_functions, g);
                const Matrix& primal_gauss_shape_derivatives =
                    primal_shape_function_gradients[g];
                typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType data;
                r_rans_primal_element.CalculateConvectionDiffusionReactionData(
                    data, primal_gauss_shape_functions,
                    primal_gauss_shape_derivatives, r_primal_process_info);
                const double scalar_value = CalculateElementScalarValue(
                    primal_gauss_shape_functions, primal_gauss_shape_derivatives,
                    r_rans_primal_element, data, r_primal_process_info);

                const double scalar_value_sensitivity =
                    (scalar_value - scalar_value_reference) / Delta;

                RansModellingApplicationTestUtilities::CheckNear(
                    adjoint_scalar_sensitivities[i_node], scalar_value_sensitivity,
                    RelativePrecision, AbsolutePrecision);

                PerturbVariable(r_node) -= Delta;
            }
        }
    }
}

template <typename TEvmElement, typename TEvmAdjointElement>
void RunScalarKEpsilon2D3NElementTest(
    const std::string PrimalElementName,
    const std::string AdjointElementName,
    const Variable<double>& rPerturbVariable,
    std::function<void(BoundedVector<double, 3>&,
                       const Vector&,
                       const Matrix&,
                       const typename TEvmAdjointElement::BaseType&,
                       const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                       const ProcessInfo&)> CalculateElementScalarValueAdjointSensitivities,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&)> CalculateElementScalarValue,
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

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.GetEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.GetEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateReactionTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateReactionTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateSourceTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateSourceTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.GetEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.GetEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateReactionTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateReactionTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateSourceTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rProcessInfo);
            noalias(rOutput) = rOutput;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            return rElement.CalculateSourceTerm(rData, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-7, 1e-6, 1e-14);
}
} // namespace Testing
} // namespace Kratos