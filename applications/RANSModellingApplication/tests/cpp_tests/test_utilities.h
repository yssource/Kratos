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
#include <cmath>
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_process.h"
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_sensitivities_process.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/checks.h"

#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
namespace RansModellingApplicationTestUtilities
{
typedef ModelPart::NodeType NodeType;
typedef ModelPart::ElementType ElementType;
typedef Geometry<NodeType> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

void IsValuesRelativelyNear(const double ValueA, const double ValueB, const double Tolerance)
{
    KRATOS_CHECK_IS_FALSE(!std::isfinite(ValueA));
    KRATOS_CHECK_IS_FALSE(!std::isfinite(ValueB));

    if (abs(ValueA) < std::numeric_limits<double>::epsilon())
    {
        KRATOS_CHECK_NEAR(ValueB, 0.0, std::numeric_limits<double>::epsilon());
        KRATOS_WARNING("IsValuesRelativelyNear")
            << "Comparing values closer to zero. ValueA/ValueB < Tolerance [ "
            << ValueA << " / " << ValueB << " < "
            << std::numeric_limits<double>::epsilon() << " ]\n";
    }
    else
    {
        const double relative_value = (1 - ValueB / ValueA);
        KRATOS_CHECK_NEAR(relative_value, 0.0, Tolerance);
    }
}

void CalculateResidual(Vector& residual, Element& rElement, ProcessInfo& rProcessInfo)
{
    const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];

    Vector rhs, nodal_scalar_values, current_nodal_scalar_rate_values,
        old_nodal_scalar_rate_values;
    Matrix damping_matrix, mass_matrix;

    rElement.CalculateRightHandSide(rhs, rProcessInfo);
    rElement.CalculateDampingMatrix(damping_matrix, rProcessInfo);
    rElement.CalculateMassMatrix(mass_matrix, rProcessInfo);

    rElement.GetFirstDerivativesVector(nodal_scalar_values);
    rElement.GetSecondDerivativesVector(current_nodal_scalar_rate_values);
    rElement.GetSecondDerivativesVector(old_nodal_scalar_rate_values, 1);

    noalias(current_nodal_scalar_rate_values) =
        current_nodal_scalar_rate_values * (1 - bossak_alpha) +
        old_nodal_scalar_rate_values * bossak_alpha;

    noalias(residual) = rhs;
    noalias(residual) -= prod(damping_matrix, nodal_scalar_values);
    noalias(residual) -= prod(mass_matrix, current_nodal_scalar_rate_values);
}

void GetElementData(Vector& rGaussWeights,
                    Matrix& rShapeFunctions,
                    ShapeFunctionDerivativesArrayType& rShapeFunctionDerivatives,
                    const ElementType& rElement)
{
    RansCalculationUtilities().CalculateGeometryData(
        rElement.GetGeometry(), rElement.GetIntegrationMethod(), rGaussWeights,
        rShapeFunctions, rShapeFunctionDerivatives);
}

void InitializeVariableWithValues(ModelPart& rModelPart,
                                  const Variable<double>& rVariable,
                                  const double Value,
                                  const std::size_t TimeStep = 0)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        r_node.FastGetSolutionStepValue(rVariable, TimeStep) = Value;
    }
}

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<double>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps)
{
    std::string seed_str = "KratosRANSModellingTestSeed_" + rVariable.Name();
    std::seed_seq seed(seed_str.begin(), seed_str.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    for (std::size_t i = 0; i < rModelPart.NumberOfNodes(); ++i)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i);
        for (std::size_t i_step = 0; i_step < TimeSteps; ++i_step)
            r_node.FastGetSolutionStepValue(rVariable, i_step) = distribution(generator);
    }
}

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps)
{
    std::string seed_str = "KratosRANSModellingTestSeed_" + rVariable.Name();
    std::seed_seq seed(seed_str.begin(), seed_str.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    for (std::size_t i = 0; i < rModelPart.NumberOfNodes(); ++i)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i);
        for (std::size_t i_step = 0; i_step < TimeSteps; ++i_step)
        {
            array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(rVariable, i_step);
            r_vector[0] = distribution(generator);
            r_vector[1] = distribution(generator);
            r_vector[2] = distribution(generator);
        }
    }
}

void RunGaussPointScalarSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    std::function<void(std::vector<double>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Vector>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance)
{
    const int number_of_elements = rModelPart.NumberOfElements();

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    rYPlusProcess.Check();

    RansCalculationUtilities rans_calculation_utilities;

    for (int i = 0; i < number_of_elements; ++i)
    {
        ElementType& r_element = *(rModelPart.ElementsBegin() + i);
        GeometryType& r_geometry = r_element.GetGeometry();
        r_element.Check(r_process_info);

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_function_derivatives;
        GetElementData(gauss_weights, shape_functions, shape_function_derivatives, r_element);

        const int number_of_gauss_points = gauss_weights.size();
        const int number_of_nodes = r_geometry.PointsNumber();

        for (int g = 0; g < number_of_gauss_points; ++g)
        {
            rYPlusProcess.Execute();
            UpdateVariablesInModelPart(rModelPart);

            const Vector& gauss_shape_functions = row(shape_functions, g);
            const Matrix& r_shape_function_derivatives = shape_function_derivatives[g];

            std::vector<Vector> analytical_sensitivities;
            CalculateSensitivities(analytical_sensitivities, r_element, gauss_shape_functions,
                                   r_shape_function_derivatives, r_process_info);

            std::vector<double> values_0;
            CalculatePrimalQuantities(values_0, r_element, gauss_shape_functions,
                                      r_shape_function_derivatives, r_process_info);

            // calculating finite difference sensitivities
            for (int i = 0; i < number_of_nodes; ++i)
            {
                NodeType& r_node = r_geometry[i];
                PerturbVariable(r_node) += Delta;

                rYPlusProcess.Execute();
                UpdateVariablesInModelPart(rModelPart);

                Vector current_gauss_weights;
                Matrix current_shape_functions;
                ShapeFunctionDerivativesArrayType current_shape_function_derivatives;
                GetElementData(current_gauss_weights, current_shape_functions,
                               current_shape_function_derivatives, r_element);

                std::vector<double> values;
                CalculatePrimalQuantities(
                    values, r_element, row(current_shape_functions, g),
                    current_shape_function_derivatives[g], r_process_info);

                for (int j = 0; j < static_cast<int>(analytical_sensitivities.size()); ++j)
                {
                    const double fd_sensitivity = ((values[j] - values_0[j]) / Delta);
                    IsValuesRelativelyNear(
                        fd_sensitivity, analytical_sensitivities[j][i], Tolerance);
                }

                PerturbVariable(r_node) -= Delta;
            }
        }
    }
}

void RunGaussPointVectorSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    std::function<void(std::vector<double>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Matrix>&, const ElementType&, const Vector&, const Matrix&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&, const int Dim)> PerturbVariable,
    const double Delta,
    const double Tolerance)
{
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    const int domain_size = r_process_info[DOMAIN_SIZE];

    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
    {
        auto calculate_sensitivities = [CalculateSensitivities, i_dim](
                                           std::vector<Vector>& riDimSensitivities,
                                           const ElementType& rElement,
                                           const Vector& rGaussShapeFunctions,
                                           const Matrix& rGaussShapeFunctionDerivatives,
                                           const ProcessInfo& rCurrentProcessInfo) {
            std::vector<Matrix> analytical_sensitivities;
            CalculateSensitivities(analytical_sensitivities, rElement, rGaussShapeFunctions,
                                   rGaussShapeFunctionDerivatives, rCurrentProcessInfo);
            for (std::size_t i_var = 0; i_var < analytical_sensitivities.size(); ++i_var)
                riDimSensitivities.push_back(column(analytical_sensitivities[i_var], i_dim));
        };

        auto perturb_variable = [PerturbVariable, i_dim](NodeType& rNode) -> double& {
            return PerturbVariable(rNode, i_dim);
        };

        RunGaussPointScalarSensitivityTest(
            rModelPart, rYPlusProcess, CalculatePrimalQuantities, calculate_sensitivities,
            UpdateVariablesInModelPart, perturb_variable, Delta, Tolerance);
    }
}

void RunElementResidualScalarSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    Process& rPrimalYPlusProcess,
    Process& rAdjointYPlusProcess,
    Process& rYPlusSensitivitiesProcess,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> CalculateElementResidualScalarSensitivity,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0)
{
    std::size_t number_of_elements = rPrimalModelPart.NumberOfElements();

    KRATOS_ERROR_IF(number_of_elements != rAdjointModelPart.NumberOfElements())
        << "Number of elements mismatch.";

    rAdjointModelPart.GetProcessInfo()[DELTA_TIME] = -1.0 * rPrimalModelPart.GetProcessInfo()[DELTA_TIME];

    // Calculate initial y_plus values
    rPrimalYPlusProcess.Check();
    rPrimalYPlusProcess.Execute();
    UpdateVariablesInModelPart(rPrimalModelPart);

    rAdjointYPlusProcess.Check();
    rAdjointYPlusProcess.Execute();
    UpdateVariablesInModelPart(rAdjointModelPart);

    // Calculate adjoint values
    rYPlusSensitivitiesProcess.Check();
    rYPlusSensitivitiesProcess.Execute();

    ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();
    ProcessInfo& r_adjoint_process_info = rAdjointModelPart.GetProcessInfo();

    const int domain_size = r_primal_process_info[DOMAIN_SIZE];
    KRATOS_ERROR_IF(domain_size != r_adjoint_process_info[DOMAIN_SIZE])
        << "Domain size mismatch.";

    Matrix adjoint_total_element_residual_sensitivity, damping_matrix, mass_matrix;

    for (std::size_t i_element = 0; i_element < number_of_elements; ++i_element)
    {
        ElementType& r_adjoint_element = *(rAdjointModelPart.ElementsBegin() + i_element);
        // r_adjoint_element.Check(r_adjoint_process_info);

        CalculateElementResidualScalarSensitivity(
            adjoint_total_element_residual_sensitivity, r_adjoint_element, r_adjoint_process_info);

        ElementType& r_primal_element = *(rPrimalModelPart.ElementsBegin() + i_element);
        GeometryType& r_primal_geometry = r_primal_element.GetGeometry();
        // r_primal_element.Check(r_primal_process_info);

        const std::size_t number_of_nodes = r_primal_geometry.PointsNumber();
        const std::size_t number_of_equations = adjoint_total_element_residual_sensitivity.size2();
        const int local_size = adjoint_total_element_residual_sensitivity.size1() / number_of_nodes;

        Vector residual(number_of_equations), residual_0(number_of_equations),
            residual_sensitivity(number_of_equations);

        CalculateResidual(residual_0, r_primal_element, r_primal_process_info);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_primal_geometry[i_node];
            PerturbVariable(r_node) += Delta;

            rPrimalYPlusProcess.Execute();
            UpdateVariablesInModelPart(rPrimalModelPart);

            CalculateResidual(residual, r_primal_element, r_primal_process_info);

            noalias(residual_sensitivity) = (residual - residual_0) / Delta;

            for (std::size_t i_check_equation = 0; i_check_equation < number_of_equations; ++i_check_equation)
            {
                const double current_adjoint_shape_sensitivity =
                    adjoint_total_element_residual_sensitivity(i_node * local_size + DerivativesOffset , i_check_equation);

                IsValuesRelativelyNear(residual_sensitivity[i_check_equation],
                                       current_adjoint_shape_sensitivity, Tolerance);
            }

            PerturbVariable(r_node) -= Delta;
        }
    }
}

void RunElementResidualVectorSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    Process& rPrimalYPlusProcess,
    Process& rAdjointYPlusProcess,
    Process& rYPlusSensitivitiesProcess,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)> CalculateElementResidualVectorSensitivity,
    std::function<double&(NodeType&, const int)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0)
{
    ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();

    const int domain_size = r_primal_process_info[DOMAIN_SIZE];

    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
    {
        auto calculate_sensitivities = [CalculateElementResidualVectorSensitivity, i_dim,
                                        domain_size, DerivativesOffset](Matrix& rDimSensitivities,
                                                     ElementType& rElement,
                                                     ProcessInfo& rCurrentProcessInfo) {
            Matrix sensitivities;
            CalculateElementResidualVectorSensitivity(sensitivities, rElement,
                                                      rCurrentProcessInfo);

            const int number_of_equations = sensitivities.size2();
            const int number_of_nodes = rElement.GetGeometry().PointsNumber();
            const int local_size = sensitivities.size1() / number_of_nodes;
            rDimSensitivities.resize(number_of_nodes, number_of_equations);

            for (int i = 0; i < number_of_nodes; ++i)
                for (int j = 0; j < number_of_equations; ++j)
                    rDimSensitivities(i, j) = sensitivities(i * local_size + i_dim + DerivativesOffset, j);
        };

        auto perturb_variable = [PerturbVariable, i_dim](NodeType& rNode) -> double& {
            return PerturbVariable(rNode, i_dim);
        };

        RunElementResidualScalarSensitivityTest(
            rPrimalModelPart, rAdjointModelPart, rPrimalYPlusProcess,
            rAdjointYPlusProcess, rYPlusSensitivitiesProcess, UpdateVariablesInModelPart,
            calculate_sensitivities, perturb_variable, Delta, Tolerance);
    }
}

void RunNodalScalarSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    std::function<void(std::vector<double>&, const NodeType&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Vector>&, const ElementType&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance)
{
    const int number_of_elements = rModelPart.NumberOfElements();

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    rYPlusProcess.Check();

    RansCalculationUtilities rans_calculation_utilities;

    for (int i = 0; i < number_of_elements; ++i)
    {
        ElementType& r_element = *(rModelPart.ElementsBegin() + i);
        GeometryType& r_geometry = r_element.GetGeometry();
        r_element.Check(r_process_info);

        const int number_of_nodes = r_geometry.PointsNumber();

        rYPlusProcess.Execute();
        UpdateVariablesInModelPart(rModelPart);

        std::vector<Vector> analytical_sensitivities;
        CalculateSensitivities(analytical_sensitivities, r_element, r_process_info);

        // calculating finite difference sensitivities
        for (int i = 0; i < number_of_nodes; ++i)
        {
            NodeType& r_node = r_geometry[i];

            std::vector<double> values_0;
            CalculatePrimalQuantities(values_0, r_node, r_process_info);

            PerturbVariable(r_node) += Delta;

            rYPlusProcess.Execute();
            UpdateVariablesInModelPart(rModelPart);

            std::vector<double> values;
            CalculatePrimalQuantities(values, r_node, r_process_info);

            for (int j = 0; j < static_cast<int>(analytical_sensitivities.size()); ++j)
            {
                const double fd_sensitivity = ((values[j] - values_0[j]) / Delta);
                IsValuesRelativelyNear(
                    fd_sensitivity, analytical_sensitivities[j][i], Tolerance);
            }

            PerturbVariable(r_node) -= Delta;
        }
    }
}

void RunNodalVectorSensitivityTest(
    ModelPart& rModelPart,
    Process& rYPlusProcess,
    std::function<void(std::vector<double>&, const NodeType&, const ProcessInfo&)> CalculatePrimalQuantities,
    std::function<void(std::vector<Matrix>&, const ElementType&, const ProcessInfo&)> CalculateSensitivities,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<double&(NodeType&, const int Dim)> PerturbVariable,
    const double Delta,
    const double Tolerance)
{
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    const int domain_size = r_process_info[DOMAIN_SIZE];

    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
    {
        auto calculate_sensitivities = [CalculateSensitivities, i_dim](
                                           std::vector<Vector>& riDimSensitivities,
                                           const ElementType& rElement,
                                           const ProcessInfo& rCurrentProcessInfo) {
            std::vector<Matrix> analytical_sensitivities;
            CalculateSensitivities(analytical_sensitivities, rElement, rCurrentProcessInfo);
            for (std::size_t i_var = 0; i_var < analytical_sensitivities.size(); ++i_var)
                riDimSensitivities.push_back(column(analytical_sensitivities[i_var], i_dim));
        };

        auto perturb_variable = [PerturbVariable, i_dim](NodeType& rNode) -> double& {
            return PerturbVariable(rNode, i_dim);
        };

        RunNodalScalarSensitivityTest(
            rModelPart, rYPlusProcess, CalculatePrimalQuantities, calculate_sensitivities,
            UpdateVariablesInModelPart, perturb_variable, Delta, Tolerance);
    }
}
} // namespace RansModellingApplicationTestUtilities
} // namespace Kratos