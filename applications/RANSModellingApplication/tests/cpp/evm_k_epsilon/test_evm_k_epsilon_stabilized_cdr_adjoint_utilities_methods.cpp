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
#include "utilities/time_discretization.h"

#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_element.h"

#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_utilities.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_utilities.h"
#include "test_evm_k_epsilon_stabilized_cdr_adjoint_interface_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateStabilizationTauScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointUtilitiesMethods)
{
    const unsigned int number_of_nodes = 3;
    const unsigned int dimension = 2;
    using adjoint_element = RansEvmKAdjoint<dimension, number_of_nodes>;
    using primal_element = RansEvmKElement<dimension, number_of_nodes>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable, dimension, number_of_nodes](
            BoundedVector<double, number_of_nodes>& rOutput,
            const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            const array_1d<double, 3>& r_velocity =
                rElement.EvaluateInPoint(VELOCITY, rShapeFunctions);
            const auto& r_parameter_derivatives =
                rElement.GetGeometryParameterDerivatives();
            BoundedMatrix<double, dimension, dimension> contravariant_metric_tensor;
            rElement.CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[GaussIndex]);

            const double reaction = rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double effective_kinematic_viscosity =
                rElement.CalculateEffectiveKinematicViscosity(
                    rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double delta_time = rElement.GetDeltaTime(rProcessInfo);

            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            const double bossak_gamma =
                TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
            const double dynamic_tau = rProcessInfo[DYNAMIC_TAU];

            double tau{0.0}, element_length{0.0};
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, r_velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            BoundedVector<double, number_of_nodes> effective_kinematic_viscosity_derivatives;
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, perturb_variable,
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);

            BoundedVector<double, number_of_nodes> reaction_derivatives;
            rElement.CalculateReactionTermScalarDerivatives(
                reaction_derivatives, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);

            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauScalarDerivatives<number_of_nodes>(
                rOutput, tau, effective_kinematic_viscosity, reaction, element_length,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            return tau;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            const array_1d<double, 3>& r_velocity =
                rElement.EvaluateInPoint(VELOCITY, rShapeFunctions);
            const auto& r_parameter_derivatives =
                rElement.GetGeometryParameterDerivatives();
            BoundedMatrix<double, dimension, dimension> contravariant_metric_tensor;
            rElement.CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[GaussIndex]);

            const double reaction = rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double effective_kinematic_viscosity =
                rElement.CalculateEffectiveKinematicViscosity(
                    rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double delta_time = rElement.GetDeltaTime(rProcessInfo);

            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            const double bossak_gamma =
                TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
            const double dynamic_tau = rProcessInfo[DYNAMIC_TAU];

            double tau{0.0}, element_length{0.0};
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, r_velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            return tau;
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateResidualScalarDerivative_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointUtilitiesMethods)
{
    const unsigned int number_of_nodes = 3;
    const unsigned int dimension = 2;
    using adjoint_element = RansEvmKAdjoint<dimension, number_of_nodes>;
    using primal_element = RansEvmKElement<dimension, number_of_nodes>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable, dimension, number_of_nodes](
            BoundedVector<double, number_of_nodes>& rOutput,
            const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            const double primal_variable_value = rElement.EvaluateInPoint(
                rElement.GetPrimalVariable(), rShapeFunctions);
            const array_1d<double, 3>& r_velocity =
                rElement.EvaluateInPoint(VELOCITY, rShapeFunctions);

            const double reaction = rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double source = rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);

            BoundedVector<double, number_of_nodes> reaction_derivatives;
            rElement.CalculateReactionTermScalarDerivatives(
                reaction_derivatives, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);

            BoundedVector<double, number_of_nodes> source_derivatives;
            rElement.CalculateSourceTermScalarDerivatives(
                source_derivatives, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);

            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualScalarDerivative(
                rOutput, primal_variable_value, reaction, r_velocity,
                reaction_derivatives, source_derivatives, rShapeFunctions,
                rShapeDerivatives, rElement.GetPrimalVariable(), perturb_variable);

            const double relaxed_scalar_rate = rElement.EvaluateInPoint(
                rElement.GetPrimalRelaxedRateVariable(), rShapeFunctions);

            array_1d<double, 3> scalar_gradient;
            rElement.CalculateGradient(
                scalar_gradient, rElement.GetPrimalVariable(), rShapeDerivatives);
            const double velocity_dot_scalar_gradient =
                inner_prod(r_velocity, scalar_gradient);

            double residual = relaxed_scalar_rate;
            residual += velocity_dot_scalar_gradient;
            residual += reaction * primal_variable_value;
            residual -= source;

            return residual;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            const Variable<double>& primal_variable = rElement.GetPrimalVariable();

            const array_1d<double, 3>& velocity =
                rElement.EvaluateInPoint(VELOCITY, rShapeFunctions);

            const double reaction = rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);

            array_1d<double, 3> variable_gradient;
            rElement.CalculateGradient(variable_gradient, primal_variable, rShapeDerivatives);
            const double velocity_dot_variable_gradient =
                inner_prod(velocity, variable_gradient);
            const double variable_value =
                rElement.EvaluateInPoint(primal_variable, rShapeFunctions);

            const double source = rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);

            const double relaxed_variable_acceleration =
                rElement.GetScalarVariableRelaxedAcceleration(rShapeFunctions);

            double residual = relaxed_variable_acceleration;
            residual += velocity_dot_variable_gradient;
            residual += reaction * variable_value;
            residual -= source;
            return residual;
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateChiScalarDerivative_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointUtilitiesMethods)
{
    const unsigned int number_of_nodes = 3;
    const unsigned int dimension = 2;
    using adjoint_element = RansEvmKAdjoint<dimension, number_of_nodes>;
    using primal_element = RansEvmKElement<dimension, number_of_nodes>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable, dimension, number_of_nodes](
            BoundedVector<double, number_of_nodes>& rOutput,
            const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            const Variable<double>& primal_variable = rElement.GetPrimalVariable();

            const array_1d<double, 3>& r_velocity =
                rElement.EvaluateInPoint(VELOCITY, rShapeFunctions);
            const auto& r_parameter_derivatives =
                rElement.GetGeometryParameterDerivatives();
            BoundedMatrix<double, dimension, dimension> contravariant_metric_tensor;
            rElement.CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[GaussIndex]);

            const double reaction = rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double effective_kinematic_viscosity =
                rElement.CalculateEffectiveKinematicViscosity(
                    rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double delta_time = rElement.GetDeltaTime(rProcessInfo);

            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            const double bossak_gamma =
                TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
            const double dynamic_tau = rProcessInfo[DYNAMIC_TAU];

            double tau{0.0}, element_length{0.0};
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, r_velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            BoundedVector<double, number_of_nodes> reaction_derivatives;
            rElement.CalculateReactionTermScalarDerivatives(
                reaction_derivatives, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);

            array_1d<double, 3> scalar_gradient;
            rElement.CalculateGradient(scalar_gradient, primal_variable, rShapeDerivatives);

            const double velocity_magnitude = norm_2(r_velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double scalar_gradient_norm = norm_2(scalar_gradient);

            double chi{0.0}, k1{0.0}, k2{0.0};

            if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);
            }

            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiScalarDerivatives(
                rOutput, chi, element_length, bossak_alpha, bossak_gamma,
                delta_time, reaction, dynamic_tau, reaction_derivatives);

            return chi;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            const Variable<double>& primal_variable = rElement.GetPrimalVariable();

            const array_1d<double, 3>& r_velocity =
                rElement.EvaluateInPoint(VELOCITY, rShapeFunctions);
            const auto& r_parameter_derivatives =
                rElement.GetGeometryParameterDerivatives();
            BoundedMatrix<double, dimension, dimension> contravariant_metric_tensor;
            rElement.CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[GaussIndex]);

            const double reaction = rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double effective_kinematic_viscosity =
                rElement.CalculateEffectiveKinematicViscosity(
                    rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            const double delta_time = rElement.GetDeltaTime(rProcessInfo);

            const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
            const double bossak_gamma =
                TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
            const double dynamic_tau = rProcessInfo[DYNAMIC_TAU];

            double tau{0.0}, element_length{0.0};
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, r_velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            array_1d<double, 3> scalar_gradient;
            rElement.CalculateGradient(scalar_gradient, primal_variable, rShapeDerivatives);

            const double velocity_magnitude = norm_2(r_velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double scalar_gradient_norm = norm_2(scalar_gradient);

            double chi{0.0}, k1{0.0}, k2{0.0};

            if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);
            }

            return chi;
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculatePositivityPreservationCoefficientScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointUtilitiesMethods)
{
    const unsigned int TNumNodes = 3;
    const unsigned int TDim = 2;
    using adjoint_element = RansEvmKAdjoint<TDim, TNumNodes>;
    using primal_element = RansEvmKElement<TDim, TNumNodes>;
    const Variable<double>& rDerivativeVariable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity = [rDerivativeVariable, TDim, TNumNodes](
                                            BoundedVector<double, TNumNodes>& rOutput,
                                            const Vector& gauss_shape_functions,
                                            const Matrix& r_shape_derivatives,
                                            const typename adjoint_element::BaseType& rElement,
                                            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& current_data,
                                            const ProcessInfo& rCurrentProcessInfo,
                                            const int g) {
        const auto& r_parameter_derivatives = rElement.GetGeometryParameterDerivatives();
        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
        array_1d<double, 3> scalar_gradient;
        const Variable<double>& primal_variable = rElement.GetPrimalVariable();

        const double delta_time = rElement.GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            scalar_convective_terms, effective_kinematic_viscosity_derivatives,
            reaction_derivatives, tau_derivatives, source_derivatives,
            s_derivatives, chi_derivatives, scalar_gradient_norm_derivative,
            residual_derivatives, absolute_residual_derivatives;

        rElement.CalculateContravariantMetricTensor(contravariant_metric_tensor,
                                                    r_parameter_derivatives[g]);

        const array_1d<double, 3>& velocity =
            rElement.EvaluateInPoint(VELOCITY, gauss_shape_functions);
        rElement.GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double scalar_value =
            rElement.EvaluateInPoint(primal_variable, gauss_shape_functions);

        rElement.CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
        rElement.GetConvectionOperator(scalar_convective_terms, scalar_gradient,
                                       r_shape_derivatives);

        const double effective_kinematic_viscosity = rElement.CalculateEffectiveKinematicViscosity(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
            effective_kinematic_viscosity_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        const double reaction = rElement.CalculateReactionTerm(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateReactionTermScalarDerivatives(
            reaction_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        const double source = rElement.CalculateSourceTerm(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateSourceTermScalarDerivatives(
            source_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        double tau, element_length;
        StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor,
            reaction, effective_kinematic_viscosity, bossak_alpha, bossak_gamma,
            delta_time, dynamic_tau);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauScalarDerivatives(
            tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
            effective_kinematic_viscosity_derivatives, reaction_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            s_derivatives, reaction, reaction_derivatives);

        const double velocity_dot_scalar_gradient = inner_prod(velocity, scalar_gradient);

        const double velocity_magnitude = norm_2(velocity);
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
        const double scalar_gradient_norm = norm_2(scalar_gradient);

        const double relaxed_scalar_rate = rElement.EvaluateInPoint(
            rElement.GetPrimalRelaxedRateVariable(), gauss_shape_functions);

        double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0},
            positivity_preserving_coeff{0.0};
        if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
            velocity_magnitude_square > std::numeric_limits<double>::epsilon())
        {
            residual = relaxed_scalar_rate;
            residual += velocity_dot_scalar_gradient;
            residual += reaction * scalar_value;
            residual -= source;
            residual = std::abs(residual);

            StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, velocity_magnitude, tau,
                effective_kinematic_viscosity, reaction, bossak_alpha,
                bossak_gamma, delta_time, element_length, dynamic_tau);

            positivity_preserving_coeff =
                residual * chi / (velocity_magnitude_square * scalar_gradient_norm);
        }

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiScalarDerivatives(
            chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
            delta_time, reaction, dynamic_tau, reaction_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarGradientScalarDerivative(
            scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualScalarDerivative(
            residual_derivatives, scalar_value, reaction, velocity,
            reaction_derivatives, source_derivatives, gauss_shape_functions,
            r_shape_derivatives, primal_variable, rDerivativeVariable);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            absolute_residual_derivatives, residual, residual_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePositivityPreservationCoefficientScalarDerivatives(
            rOutput, chi, residual, scalar_gradient_norm, velocity_magnitude_square,
            chi_derivatives, absolute_residual_derivatives,
            scalar_gradient_norm_derivative, primal_variable, rDerivativeVariable);

        return positivity_preserving_coeff;
    };

    auto calculate_scalar_value = [TDim, TNumNodes](
                                      const Vector& gauss_shape_functions,
                                      const Matrix& r_shape_derivatives,
                                      const typename primal_element::BaseType& rElement,
                                      const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& r_current_data,
                                      const ProcessInfo& rCurrentProcessInfo, const int g) {
        const auto& r_parameter_derivatives = rElement.GetGeometryParameterDerivatives();
        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
        array_1d<double, 3> variable_gradient;
        const Variable<double>& primal_variable = rElement.GetPrimalVariable();

        const double delta_time = rElement.GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        rElement.CalculateContravariantMetricTensor(contravariant_metric_tensor,
                                                    r_parameter_derivatives[g]);

        const array_1d<double, 3>& velocity =
            rElement.EvaluateInPoint(VELOCITY, gauss_shape_functions);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        rElement.GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
        const double velocity_magnitude = norm_2(velocity);

        const double effective_kinematic_viscosity = rElement.CalculateEffectiveKinematicViscosity(
            r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        const double variable_gradient_norm =
            rElement.GetScalarVariableGradientNorm(r_shape_derivatives);
        const double relaxed_variable_acceleration =
            rElement.GetScalarVariableRelaxedAcceleration(gauss_shape_functions);
        rElement.CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);

        const double reaction = rElement.CalculateReactionTerm(
            r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        double tau, element_length;
        StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor,
            reaction, effective_kinematic_viscosity, bossak_alpha, bossak_gamma,
            delta_time, dynamic_tau);

        // Calculate residual for cross wind dissipation coefficient
        double positivity_preserving_coefficient{0.0}, k1{0.0}, k2{0.0}, chi{0.0};
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

        const double velocity_dot_variable_gradient =
            inner_prod(velocity, variable_gradient);
        const double variable_value =
            rElement.EvaluateInPoint(primal_variable, gauss_shape_functions);

        if (variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
            velocity_magnitude_square > std::numeric_limits<double>::epsilon())
        {
            const double source = rElement.CalculateSourceTerm(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double residual = relaxed_variable_acceleration;
            residual += velocity_dot_variable_gradient;
            residual += reaction * variable_value;
            residual -= source;
            residual = std::abs(residual);

            StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, velocity_magnitude, tau,
                effective_kinematic_viscosity, reaction, bossak_alpha,
                bossak_gamma, delta_time, element_length, dynamic_tau);

            positivity_preserving_coefficient =
                residual * chi / (velocity_magnitude_square * variable_gradient_norm);
        }

        return positivity_preserving_coefficient;
    };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", rDerivativeVariable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateCrossWindDiffusionCoeffScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointUtilitiesMethods)
{
    const unsigned int TNumNodes = 3;
    const unsigned int TDim = 2;
    using adjoint_element = RansEvmKAdjoint<TDim, TNumNodes>;
    using primal_element = RansEvmKElement<TDim, TNumNodes>;
    const Variable<double>& rDerivativeVariable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity = [rDerivativeVariable, TDim, TNumNodes](
                                            BoundedVector<double, TNumNodes>& rOutput,
                                            const Vector& gauss_shape_functions,
                                            const Matrix& r_shape_derivatives,
                                            const typename adjoint_element::BaseType& rElement,
                                            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& current_data,
                                            const ProcessInfo& rCurrentProcessInfo,
                                            const int g) {
        const auto& r_parameter_derivatives = rElement.GetGeometryParameterDerivatives();
        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
        array_1d<double, 3> scalar_gradient;
        const Variable<double>& primal_variable = rElement.GetPrimalVariable();

        const double delta_time = rElement.GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            scalar_convective_terms, effective_kinematic_viscosity_derivatives,
            reaction_derivatives, tau_derivatives, source_derivatives, s_derivatives,
            chi_derivatives, scalar_gradient_norm_derivative, residual_derivatives,
            absolute_residual_derivatives, positivity_preserving_coeff_derivatives,
            absolute_reaction_tilde_derivatives, psi_one_derivatives, psi_two_derivatives;

        rElement.CalculateContravariantMetricTensor(contravariant_metric_tensor,
                                                    r_parameter_derivatives[g]);

        const array_1d<double, 3>& velocity =
            rElement.EvaluateInPoint(VELOCITY, gauss_shape_functions);
        rElement.GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double scalar_value =
            rElement.EvaluateInPoint(primal_variable, gauss_shape_functions);

        rElement.CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
        rElement.GetConvectionOperator(scalar_convective_terms, scalar_gradient,
                                       r_shape_derivatives);

        const double effective_kinematic_viscosity = rElement.CalculateEffectiveKinematicViscosity(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
            effective_kinematic_viscosity_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        const double reaction = rElement.CalculateReactionTerm(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateReactionTermScalarDerivatives(
            reaction_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        const double source = rElement.CalculateSourceTerm(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateSourceTermScalarDerivatives(
            source_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        double tau, element_length;
        StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor,
            reaction, effective_kinematic_viscosity, bossak_alpha, bossak_gamma,
            delta_time, dynamic_tau);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauScalarDerivatives(
            tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
            effective_kinematic_viscosity_derivatives, reaction_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            s_derivatives, reaction, reaction_derivatives);

        const double velocity_dot_scalar_gradient = inner_prod(velocity, scalar_gradient);

        const double velocity_magnitude = norm_2(velocity);
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
        const double scalar_gradient_norm = norm_2(scalar_gradient);

        const double relaxed_scalar_rate = rElement.EvaluateInPoint(
            rElement.GetPrimalRelaxedRateVariable(), gauss_shape_functions);

        double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0};
        if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
            velocity_magnitude_square > std::numeric_limits<double>::epsilon())
        {
            residual = relaxed_scalar_rate;
            residual += velocity_dot_scalar_gradient;
            residual += reaction * scalar_value;
            residual -= source;
            residual = std::abs(residual);

            StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, velocity_magnitude, tau,
                effective_kinematic_viscosity, reaction, bossak_alpha,
                bossak_gamma, delta_time, element_length, dynamic_tau);
        }

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiScalarDerivatives(
            chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
            delta_time, reaction, dynamic_tau, reaction_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarGradientScalarDerivative(
            scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualScalarDerivative(
            residual_derivatives, scalar_value, reaction, velocity,
            reaction_derivatives, source_derivatives, gauss_shape_functions,
            r_shape_derivatives, primal_variable, rDerivativeVariable);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            absolute_residual_derivatives, residual, residual_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePositivityPreservationCoefficientScalarDerivatives(
            positivity_preserving_coeff_derivatives, chi, residual, scalar_gradient_norm,
            velocity_magnitude_square, chi_derivatives, absolute_residual_derivatives,
            scalar_gradient_norm_derivative, primal_variable, rDerivativeVariable);

        const double reaction_tilde =
            reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

        const double psi_one = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
            velocity_magnitude, tau, reaction_tilde);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiOneScalarDerivatives(
            psi_one_derivatives, velocity_magnitude, reaction_tilde, tau,
            tau_derivatives, absolute_reaction_tilde_derivatives);

        const double psi_two = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
            reaction_tilde, tau, element_length);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiTwoScalarDerivatives(
            psi_two_derivatives, element_length, tau, reaction_tilde, tau_derivatives,
            reaction_derivatives, absolute_reaction_tilde_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateCrossWindDiffusionCoeffScalarDerivatives(
            rOutput, k2, psi_one, psi_two, element_length,
            effective_kinematic_viscosity, psi_one_derivatives,
            psi_two_derivatives, effective_kinematic_viscosity_derivatives);

        return k2;
    };

    auto calculate_scalar_value =
        [TDim, TNumNodes](
            const Vector& gauss_shape_functions, const Matrix& r_shape_derivatives,
            const typename primal_element::BaseType& rElement,
            const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& r_current_data,
            const ProcessInfo& rCurrentProcessInfo, const int g) {
            const auto& r_parameter_derivatives =
                rElement.GetGeometryParameterDerivatives();
            BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
            array_1d<double, 3> variable_gradient;
            const Variable<double>& primal_variable = rElement.GetPrimalVariable();

            const double delta_time = rElement.GetDeltaTime(rCurrentProcessInfo);
            const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
            const double bossak_gamma =
                TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
            const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

            rElement.CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const array_1d<double, 3>& velocity =
                rElement.EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            rElement.GetConvectionOperator(velocity_convective_terms, velocity,
                                           r_shape_derivatives);
            const double velocity_magnitude = norm_2(velocity);

            const double effective_kinematic_viscosity = rElement.CalculateEffectiveKinematicViscosity(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            const double variable_gradient_norm =
                rElement.GetScalarVariableGradientNorm(r_shape_derivatives);
            const double relaxed_variable_acceleration =
                rElement.GetScalarVariableRelaxedAcceleration(gauss_shape_functions);
            rElement.CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);

            const double reaction = rElement.CalculateReactionTerm(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            // Calculate residual for cross wind dissipation coefficient
            double k1{0.0}, k2{0.0}, chi{0.0};
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

            const double velocity_dot_variable_gradient =
                inner_prod(velocity, variable_gradient);
            const double variable_value =
                rElement.EvaluateInPoint(primal_variable, gauss_shape_functions);

            if (variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                const double source = rElement.CalculateSourceTerm(
                    r_current_data, gauss_shape_functions, r_shape_derivatives,
                    rCurrentProcessInfo);

                double residual = relaxed_variable_acceleration;
                residual += velocity_dot_variable_gradient;
                residual += reaction * variable_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);
            }

            return k2;
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", rDerivativeVariable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateStreamLineDiffusionCoeffScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointUtilitiesMethods)
{
    const unsigned int TNumNodes = 3;
    const unsigned int TDim = 2;
    using adjoint_element = RansEvmKAdjoint<TDim, TNumNodes>;
    using primal_element = RansEvmKElement<TDim, TNumNodes>;
    const Variable<double>& rDerivativeVariable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity = [rDerivativeVariable, TDim, TNumNodes](
                                            BoundedVector<double, TNumNodes>& rOutput,
                                            const Vector& gauss_shape_functions,
                                            const Matrix& r_shape_derivatives,
                                            const typename adjoint_element::BaseType& rElement,
                                            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& current_data,
                                            const ProcessInfo& rCurrentProcessInfo,
                                            const int g) {
        const auto& r_parameter_derivatives = rElement.GetGeometryParameterDerivatives();
        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
        array_1d<double, 3> scalar_gradient;
        const Variable<double>& primal_variable = rElement.GetPrimalVariable();

        const double delta_time = rElement.GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            scalar_convective_terms, effective_kinematic_viscosity_derivatives,
            reaction_derivatives, tau_derivatives, source_derivatives, s_derivatives,
            chi_derivatives, scalar_gradient_norm_derivative, residual_derivatives,
            absolute_residual_derivatives, positivity_preserving_coeff_derivatives,
            absolute_reaction_tilde_derivatives, psi_one_derivatives, psi_two_derivatives;

        rElement.CalculateContravariantMetricTensor(contravariant_metric_tensor,
                                                    r_parameter_derivatives[g]);

        const array_1d<double, 3>& velocity =
            rElement.EvaluateInPoint(VELOCITY, gauss_shape_functions);
        rElement.GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double scalar_value =
            rElement.EvaluateInPoint(primal_variable, gauss_shape_functions);

        rElement.CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
        rElement.GetConvectionOperator(scalar_convective_terms, scalar_gradient,
                                       r_shape_derivatives);

        const double effective_kinematic_viscosity = rElement.CalculateEffectiveKinematicViscosity(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
            effective_kinematic_viscosity_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        const double reaction = rElement.CalculateReactionTerm(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateReactionTermScalarDerivatives(
            reaction_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        const double source = rElement.CalculateSourceTerm(
            current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
        rElement.CalculateSourceTermScalarDerivatives(
            source_derivatives, rDerivativeVariable, current_data,
            gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

        double tau, element_length;
        StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor,
            reaction, effective_kinematic_viscosity, bossak_alpha, bossak_gamma,
            delta_time, dynamic_tau);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauScalarDerivatives(
            tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
            effective_kinematic_viscosity_derivatives, reaction_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            s_derivatives, reaction, reaction_derivatives);

        const double velocity_dot_scalar_gradient = inner_prod(velocity, scalar_gradient);

        const double velocity_magnitude = norm_2(velocity);
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
        const double scalar_gradient_norm = norm_2(scalar_gradient);

        const double relaxed_scalar_rate = rElement.EvaluateInPoint(
            rElement.GetPrimalRelaxedRateVariable(), gauss_shape_functions);

        double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0};
        if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
            velocity_magnitude_square > std::numeric_limits<double>::epsilon())
        {
            residual = relaxed_scalar_rate;
            residual += velocity_dot_scalar_gradient;
            residual += reaction * scalar_value;
            residual -= source;
            residual = std::abs(residual);

            StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, velocity_magnitude, tau,
                effective_kinematic_viscosity, reaction, bossak_alpha,
                bossak_gamma, delta_time, element_length, dynamic_tau);
        }

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiScalarDerivatives(
            chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
            delta_time, reaction, dynamic_tau, reaction_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarGradientScalarDerivative(
            scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualScalarDerivative(
            residual_derivatives, scalar_value, reaction, velocity,
            reaction_derivatives, source_derivatives, gauss_shape_functions,
            r_shape_derivatives, primal_variable, rDerivativeVariable);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            absolute_residual_derivatives, residual, residual_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePositivityPreservationCoefficientScalarDerivatives(
            positivity_preserving_coeff_derivatives, chi, residual, scalar_gradient_norm,
            velocity_magnitude_square, chi_derivatives, absolute_residual_derivatives,
            scalar_gradient_norm_derivative, primal_variable, rDerivativeVariable);

        const double reaction_tilde =
            reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
            absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

        const double psi_one = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
            velocity_magnitude, tau, reaction_tilde);
        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiOneScalarDerivatives(
            psi_one_derivatives, velocity_magnitude, reaction_tilde, tau,
            tau_derivatives, absolute_reaction_tilde_derivatives);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiTwoScalarDerivatives(
            psi_two_derivatives, element_length, tau, reaction_tilde, tau_derivatives,
            reaction_derivatives, absolute_reaction_tilde_derivatives);

        const double psi_two = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
            reaction_tilde, tau, element_length);

        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStreamLineDiffusionCoeffScalarDerivatives(
            rOutput, k1, element_length, tau, velocity_magnitude, reaction_tilde,
            psi_one, psi_two, psi_one_derivatives, psi_two_derivatives, tau_derivatives,
            reaction_derivatives, effective_kinematic_viscosity_derivatives);

        return k1;
    };

    auto calculate_scalar_value =
        [TDim, TNumNodes](
            const Vector& gauss_shape_functions, const Matrix& r_shape_derivatives,
            const typename primal_element::BaseType& rElement,
            const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& r_current_data,
            const ProcessInfo& rCurrentProcessInfo, const int g) {
            const auto& r_parameter_derivatives =
                rElement.GetGeometryParameterDerivatives();
            BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
            array_1d<double, 3> variable_gradient;
            const Variable<double>& primal_variable = rElement.GetPrimalVariable();

            const double delta_time = rElement.GetDeltaTime(rCurrentProcessInfo);
            const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
            const double bossak_gamma =
                TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
            const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

            rElement.CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const array_1d<double, 3>& velocity =
                rElement.EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            rElement.GetConvectionOperator(velocity_convective_terms, velocity,
                                           r_shape_derivatives);
            const double velocity_magnitude = norm_2(velocity);

            const double effective_kinematic_viscosity = rElement.CalculateEffectiveKinematicViscosity(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            const double variable_gradient_norm =
                rElement.GetScalarVariableGradientNorm(r_shape_derivatives);
            const double relaxed_variable_acceleration =
                rElement.GetScalarVariableRelaxedAcceleration(gauss_shape_functions);
            rElement.CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);

            const double reaction = rElement.CalculateReactionTerm(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            // Calculate residual for cross wind dissipation coefficient
            double k1{0.0}, k2{0.0}, chi{0.0};
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

            const double velocity_dot_variable_gradient =
                inner_prod(velocity, variable_gradient);
            const double variable_value =
                rElement.EvaluateInPoint(primal_variable, gauss_shape_functions);

            if (variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                const double source = rElement.CalculateSourceTerm(
                    r_current_data, gauss_shape_functions, r_shape_derivatives,
                    rCurrentProcessInfo);

                double residual = relaxed_variable_acceleration;
                residual += velocity_dot_variable_gradient;
                residual += reaction * variable_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);
            }

            return k1;
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", rDerivativeVariable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

} // namespace Testing

} // namespace Kratos
