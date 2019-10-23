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

#include "custom_elements/evm_k_epsilon/rans_evm_k_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_element.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_utilities.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_utilities.h"
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
                       const ProcessInfo&)> CalculateElementAdjointSensitivities)
{
    const int number_of_primal_elements = rPrimalModelPart.NumberOfElements();
    const int number_of_adjoint_elements = rAdjointModelPart.NumberOfElements();

    KRATOS_ERROR_IF(number_of_primal_elements != number_of_adjoint_elements)
        << "Primal and Adjoint model part number of elements mismatch.\n";

    for (Process* p_process : rAdjointProcessList)
        p_process->Check();

    for (Process* p_process : rAdjointProcessList)
        p_process->ExecuteInitialize();

    for (Process* p_process : rAdjointProcessList)
        p_process->ExecuteInitializeSolutionStep();

    const ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();
    const ProcessInfo& r_adjoint_process_info = rAdjointModelPart.GetProcessInfo();

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

        auto& r_adjoint_geometry = r_adjoint_element.GetGeometry();
        const Matrix& adjoint_shape_functions =
            r_adjoint_geometry.ShapeFunctionsValues(r_adjoint_integration_method);
        typename Element::GeometryType::ShapeFunctionsGradientsType adjoint_shape_function_gradients;
        r_adjoint_geometry.ShapeFunctionsIntegrationPointsGradients(
            adjoint_shape_function_gradients, r_adjoint_integration_method);

        const int number_of_gauss_points = adjoint_shape_functions.size1();
        for (int g = 0; g < number_of_gauss_points; ++g)
        {
            const Vector& r_adjoint_gauss_shape_functions =
                row(adjoint_shape_functions, g);
            const Matrix& r_adjoint_gauss_shape_derivatives =
                adjoint_shape_function_gradients[g];

            typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType adjoint_data;
            auto& r_rans_adjoint_element =
                dynamic_cast<const typename TEvmAdjointElement::BaseType&>(r_adjoint_element);
            r_rans_adjoint_element.CalculateElementData(
                adjoint_data, r_adjoint_gauss_shape_functions,
                r_adjoint_gauss_shape_derivatives, r_adjoint_process_info);

            BoundedVector<double, TNumNodes> adjoint_scalar_sensitivities;
            CalculateElementAdjointSensitivities(
                adjoint_scalar_sensitivities, r_adjoint_gauss_shape_functions,
                r_adjoint_gauss_shape_derivatives, r_rans_adjoint_element,
                adjoint_data, r_adjoint_process_info);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmKAdjoint2D3N");

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

    auto perturbation = [perturb_variable](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(perturb_variable);
    };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo) {
            //    return rElement.GetEffectiveKinematicViscosity(rData, )
        };

    RunScalarSensitivityTest<primal_element, adjoint_element, 2, 3>(
        r_primal_model_part, r_adjoint_model_part, primal_processes_list,
        adjoint_processes_list, calculate_scalar_sensitivity);
}
} // namespace Testing
} // namespace Kratos