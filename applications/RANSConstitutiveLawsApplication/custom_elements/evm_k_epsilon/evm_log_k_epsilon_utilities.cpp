#include "evm_log_k_epsilon_utilities.h"
#include <cmath>
#include <iostream>
#include <limits>

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace EvmLogKepsilonModelUtilities
{
double CalculateTurbulentViscosity(const double C_mu,
                                   const double log_turbulent_kinetic_energy,
                                   const double log_turbulent_energy_dissipation_rate,
                                   const double f_mu)
{
    return C_mu * f_mu * std::exp(2.0 * log_turbulent_kinetic_energy - log_turbulent_energy_dissipation_rate);
}

double CalculateFmu(const double y_plus)
{
    return 1.0 - std::exp(-0.0115 * y_plus);
}

double CalculateF2(const double log_turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double log_turbulent_energy_dissipation_rate)
{
    const double Re_t = std::exp(2.0 * log_turbulent_kinetic_energy -
                                 log_turbulent_energy_dissipation_rate) /
                        kinematic_viscosity;
    const double f2 = 1.0 - 0.22 * std::exp(-1.0 * std::pow(Re_t / 6.0, 2));

    return f2;
}

void CalculateStabilizationTau(double& tau,
                               double& element_length,
                               const array_1d<double, 3>& rVelocity,
                               const Matrix& rContravariantMetricTensor,
                               const double effective_kinematic_viscosity,
                               const double delta_time)
{
    unsigned int dim = rContravariantMetricTensor.size2();
    Vector velocity(dim);
    for (unsigned int d = 0; d < dim; ++d)
        velocity[d] = rVelocity[d];
    Vector temp(dim);
    noalias(temp) = prod(rContravariantMetricTensor, velocity);
    const double stab_convection = inner_prod(velocity, temp);
    const double stab_diffusion = std::pow(3.0 * effective_kinematic_viscosity, 2) *
                                  norm_frobenius(rContravariantMetricTensor);
    const double stab_dynamics = std::pow(2.0 / delta_time, 2);
    PrintIfVariableIsNegative(stab_convection);

    const double velocity_norm = norm_2(rVelocity);
    tau = 1.0 / std::sqrt(stab_dynamics + stab_convection + stab_diffusion);
    element_length = 2.0 * velocity_norm / std::sqrt(stab_convection);
}

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity)
{
    const double velocity_divergence =
        CalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;

    noalias(reynolds_stress_tensor) =
        turbulent_kinematic_viscosity *
        (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity);

    double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            source += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    return source;
}

double CalculateGamma(const double log_turbulent_kinetic_energy,
                      const double log_turbulent_energy_dissipation_rate)
{
    return std::exp(log_turbulent_energy_dissipation_rate - log_turbulent_kinetic_energy);
}

void CalculateCrossWindDiffusionParameters(double& chi,
                                           double& k1,
                                           double& k2,
                                           const double velocity_magnitude,
                                           const double tau,
                                           const double effective_kinematic_viscosity,
                                           const double reaction,
                                           const double element_length)
{
    chi = 2.0 / (std::abs(reaction) * element_length + 2.0 * velocity_magnitude);

    double value = 0.0;
    value = std::abs(0.5 * (velocity_magnitude - tau * velocity_magnitude * reaction +
                            tau * velocity_magnitude * std::abs(reaction))) *
                element_length -
            (effective_kinematic_viscosity + tau * std::pow(velocity_magnitude, 2)) +
            (reaction + tau * reaction * std::abs(reaction)) *
                std::pow(element_length, 2) / 6.0;
    k1 = std::max(value, 0.0);

    value = std::abs(0.5 * (velocity_magnitude + tau * velocity_magnitude * std::abs(reaction))) *
                element_length -
            effective_kinematic_viscosity +
            (reaction + tau * reaction * std::abs(reaction)) *
                std::pow(element_length, 2) / 6.0;
    k2 = std::max(value, 0.0);
}

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double);

} // namespace EvmLogKepsilonModelUtilities

///@}

} // namespace Kratos
