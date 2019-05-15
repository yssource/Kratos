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

namespace Kratos
{
namespace RansModellingApplicationTestUtilities
{
typedef ModelPart::NodeType NodeType;
typedef ModelPart::ElementType ElementType;
typedef Geometry<NodeType> GeometryType;

void IsValidValue(const double Value, const double Tolerance)
{
    KRATOS_CHECK_IS_FALSE(!std::isfinite(Value));
    KRATOS_CHECK_IS_FALSE(std::abs(Value) < Tolerance);
}

void IsValuesRelativelyNear(const double ValueA, const double ValueB, const double Tolerance)
{
    IsValidValue(ValueA, Tolerance);
    IsValidValue(ValueB, Tolerance);

    const double relative_value = (1 - ValueB / ValueA);

    KRATOS_CHECK_NEAR(relative_value, 0.0, Tolerance);
}

void CalculateResidual(Vector& residual,
                       Element& rElement,
                       ProcessInfo& rProcessInfo,
                       const Variable<double>& rScalarVariable,
                       const Variable<double>& rRelaxedScalarRateVariable)
{
    RansVariableUtils rans_variable_utils;

    Vector rhs, nodal_scalar_values, nodal_scalar_relaxed_rate_values;
    Matrix damping_matrix, mass_matrix;

    rElement.CalculateRightHandSide(rhs, rProcessInfo);
    rElement.CalculateDampingMatrix(damping_matrix, rProcessInfo);
    rElement.CalculateMassMatrix(mass_matrix, rProcessInfo);

    rans_variable_utils.GetNodalArray(nodal_scalar_values, rElement, rScalarVariable);
    rans_variable_utils.GetNodalArray(nodal_scalar_relaxed_rate_values,
                                      rElement, rRelaxedScalarRateVariable);

    noalias(residual) = rhs;
    noalias(residual) -= prod(damping_matrix, nodal_scalar_values);
    noalias(residual) -= prod(mass_matrix, nodal_scalar_relaxed_rate_values);
}

void RunResidualVectorSensitivityTest(
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rPrimalRelaxedRateVariable,
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    Process& rPrimalYPlusProcess,
    Process& rAdjointYPlusProcess,
    Process& rYPlusSensitivitiesProcess,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, Element&, const ProcessInfo&)> CalculateSensitivityMatrix,
    std::function<void(NodeType&, const int, const double)> PerturbVariable,
    const double Delta,
    const double Tolerance)
{
    std::size_t number_of_elements = rPrimalModelPart.NumberOfElements();

    KRATOS_ERROR_IF(number_of_elements != rAdjointModelPart.NumberOfElements())
        << "Number of elements mismatch.";

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
        CalculateSensitivityMatrix(adjoint_total_element_residual_sensitivity,
                                   r_adjoint_element, r_adjoint_process_info);

        ElementType& r_primal_element = *(rPrimalModelPart.ElementsBegin() + i_element);
        GeometryType& r_primal_geometry = r_primal_element.GetGeometry();

        const std::size_t number_of_nodes = r_primal_geometry.PointsNumber();

        Vector residual(number_of_nodes), residual_0(number_of_nodes),
            residual_sensitivity(number_of_nodes);

        CalculateResidual(residual_0, r_primal_element, r_primal_process_info,
                          rPrimalVariable, rPrimalRelaxedRateVariable);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_primal_geometry[i_node];
            for (int i_dim = 0; i_dim < domain_size; ++i_dim)
            {
                PerturbVariable(r_node, i_dim, Delta);

                rPrimalYPlusProcess.Execute();
                UpdateVariablesInModelPart(rPrimalModelPart);

                CalculateResidual(residual, r_primal_element, r_primal_process_info,
                                  rPrimalVariable, rPrimalRelaxedRateVariable);

                noalias(residual_sensitivity) = (residual - residual_0) / Delta;

                for (std::size_t i_check_node = 0; i_check_node < number_of_nodes; ++i_check_node)
                {
                    const double current_adjoint_shape_sensitivity =
                        adjoint_total_element_residual_sensitivity(
                            i_node * 2 + i_dim, i_check_node);

                    IsValuesRelativelyNear(residual_sensitivity[i_check_node],
                                           current_adjoint_shape_sensitivity, Tolerance);
                }

                PerturbVariable(r_node, i_dim, -Delta);
            }
        }
    }
}

void RunResidualScalarSensitivityTest(
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rPrimalRelaxedRateVariable,
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    Process& rPrimalYPlusProcess,
    Process& rAdjointYPlusProcess,
    Process& rYPlusSensitivitiesProcess,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, Element&, ProcessInfo&)> CalculateSensitivityMatrix,
    std::function<void(NodeType&, const double)> PerturbVariable,
    const double Delta,
    const double Tolerance)
{
    std::size_t number_of_elements = rPrimalModelPart.NumberOfElements();

    KRATOS_ERROR_IF(number_of_elements != rAdjointModelPart.NumberOfElements())
        << "Number of elements mismatch.";

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
        CalculateSensitivityMatrix(adjoint_total_element_residual_sensitivity,
                                   r_adjoint_element, r_adjoint_process_info);

        ElementType& r_primal_element = *(rPrimalModelPart.ElementsBegin() + i_element);
        GeometryType& r_primal_geometry = r_primal_element.GetGeometry();

        const std::size_t number_of_nodes = r_primal_geometry.PointsNumber();

        Vector residual(number_of_nodes), residual_0(number_of_nodes),
            residual_sensitivity(number_of_nodes);

        CalculateResidual(residual_0, r_primal_element, r_primal_process_info,
                          rPrimalVariable, rPrimalRelaxedRateVariable);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_primal_geometry[i_node];
            PerturbVariable(r_node, Delta);

            rPrimalYPlusProcess.Execute();
            UpdateVariablesInModelPart(rPrimalModelPart);

            CalculateResidual(residual, r_primal_element, r_primal_process_info,
                              rPrimalVariable, rPrimalRelaxedRateVariable);

            noalias(residual_sensitivity) = (residual - residual_0) / Delta;

            for (std::size_t i_check_node = 0; i_check_node < number_of_nodes; ++i_check_node)
            {
                const double current_adjoint_shape_sensitivity =
                    adjoint_total_element_residual_sensitivity(i_node, i_check_node);

                IsValuesRelativelyNear(residual_sensitivity[i_check_node],
                                       current_adjoint_shape_sensitivity, Tolerance);
            }

            PerturbVariable(r_node, -Delta);
        }
    }
}
} // namespace RansModellingApplicationTestUtilities
} // namespace Kratos