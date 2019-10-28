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

#include "custom_utilities/test_utilities.h"

#if !defined(KRATOS_RANS_TEST_STABILIZED_CDR_ADJOINT_INTERFACE_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_TEST_STABILIZED_CDR_ADJOINT_INTERFACE_UTILITIES_H_INCLUDED

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
    std::function<double(BoundedVector<double, TNumNodes>&,
                         const Vector&,
                         const Matrix&,
                         const typename TEvmAdjointElement::BaseType&,
                         const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValueAdjointScalarSensitivities,
    std::function<double&(ModelPart::NodeType&)> PerturbVariable,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValue,
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

    const int primal_domain_size = r_primal_process_info[DOMAIN_SIZE];
    const int adjoint_domain_size = r_adjoint_process_info[DOMAIN_SIZE];

    KRATOS_ERROR_IF(primal_domain_size != TDim)
        << "Primal model part domain size mismatch.";
    KRATOS_ERROR_IF(adjoint_domain_size != TDim)
        << "Adjoint model part domain size mismatch.";

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
            const double adjoint_scalar_value =
                CalculateElementScalarValueAdjointScalarSensitivities(
                    adjoint_scalar_sensitivities, r_adjoint_gauss_shape_functions,
                    r_adjoint_gauss_shape_derivatives, r_rans_adjoint_element,
                    adjoint_data, r_adjoint_process_info, g);

            // KRATOS_WATCH(adjoint_scalar_sensitivities);

            // Calculating reference scalar value
            execute_primal_process();
            UpdateVariablesInModelPart(rPrimalModelPart);
            const Matrix& primal_shape_functions = calculate_shape_functions(
                primal_shape_function_gradients, r_primal_element);
            const Vector& primal_gauss_shape_functions = row(primal_shape_functions, g);
            const Matrix& primal_gauss_shape_derivatives =
                primal_shape_function_gradients[g];
            typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType data;
            r_rans_primal_element.CalculateElementData(
                data, primal_gauss_shape_functions,
                primal_gauss_shape_derivatives, r_primal_process_info);

            // Check scalar value reference calculation
            const double scalar_value_reference = CalculateElementScalarValue(
                primal_gauss_shape_functions, primal_gauss_shape_derivatives,
                r_rans_primal_element, data, r_primal_process_info, g);

            KRATOS_CHECK_NEAR(adjoint_scalar_value, scalar_value_reference, 1e-12);

            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            {
                ModelPart::NodeType& r_node = r_primal_geometry[i_node];
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
                r_rans_primal_element.CalculateElementData(
                    data, primal_gauss_shape_functions,
                    primal_gauss_shape_derivatives, r_primal_process_info);
                const double scalar_value = CalculateElementScalarValue(
                    primal_gauss_shape_functions, primal_gauss_shape_derivatives,
                    r_rans_primal_element, data, r_primal_process_info, g);

                const double scalar_value_sensitivity =
                    (scalar_value - scalar_value_reference) / Delta;

                RansModellingApplicationTestUtilities::CheckNear(
                    adjoint_scalar_sensitivities[i_node], scalar_value_sensitivity,
                    RelativePrecision, AbsolutePrecision);

                // KRATOS_WATCH(scalar_value_sensitivity);

                PerturbVariable(r_node) -= Delta;
            }
        }
    }
}

template <typename TEvmElement, typename TEvmAdjointElement, unsigned int TDim, unsigned int TNumNodes>
void RunVectorSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcessList,
    std::vector<Process*>& rAdjointProcessList,
    std::function<double(BoundedMatrix<double, TNumNodes, TDim>&,
                         const Vector&,
                         const Matrix&,
                         const typename TEvmAdjointElement::BaseType&,
                         const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValueAdjointVectorSensitivities,
    std::function<double&(ModelPart::NodeType&, const int)> PerturbVariable,
    std::function<double(const Vector&,
                         const Matrix&,
                         const typename TEvmElement::BaseType&,
                         const typename TEvmElement::BaseType::ConvectionDiffusionReactionDataType&,
                         const ProcessInfo&,
                         const int GaussIndex)> CalculateElementScalarValue,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    const double Delta,
    const double RelativePrecision,
    const double AbsolutePrecision)
{
    for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
    {
        auto calculate_sensitivities =
            [CalculateElementScalarValueAdjointVectorSensitivities, i_dim](
                BoundedVector<double, TNumNodes>& rOutput,
                const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives,
                const typename TEvmAdjointElement::BaseType& rElement,
                const typename TEvmAdjointElement::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
                const ProcessInfo& rCurrentProcessInfo, const int GaussIndex) {
                BoundedMatrix<double, TNumNodes, TDim> scalar_vector_sensitivities;
                const double scalar_value = CalculateElementScalarValueAdjointVectorSensitivities(
                    scalar_vector_sensitivities, rShapeFunctions, rShapeFunctionDerivatives,
                    rElement, rData, rCurrentProcessInfo, GaussIndex);

                noalias(rOutput) = column(scalar_vector_sensitivities, i_dim);
                return scalar_value;
            };

        auto perturb_variable = [PerturbVariable,
                                 i_dim](ModelPart::NodeType& rNode) -> double& {
            return PerturbVariable(rNode, i_dim);
        };

        // RunScalarSensitivityTest<TEvmElement, TEvmAdjointElement, TDim, TNumNodes>(
        //     rPrimalModelPart, rAdjointModelPart, rPrimalProcessList, rAdjointProcessList,
        //     calculate_sensitivities, perturb_variable, CalculateElementScalarValue,
        //     UpdateVariablesInModelPart, Delta, RelativePrecision, AbsolutePrecision);
    }
}
} // namespace Testing
} // namespace Kratos
#endif