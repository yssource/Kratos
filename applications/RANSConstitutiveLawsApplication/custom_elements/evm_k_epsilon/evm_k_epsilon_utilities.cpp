#include "evm_k_epsilon_utilities.h"
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

namespace EvmKepsilonModelUtilities
{
double CalculateTurbulentViscosity(const double C_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double f_mu)
{
    return C_mu * f_mu * std::pow(turbulent_kinetic_energy, 2) / turbulent_energy_dissipation_rate;
}

double CalculateFmu(const double y_plus)
{
    return 1.0 - std::exp(-0.0115 * y_plus);
}

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate)
{
    if (turbulent_energy_dissipation_rate == 0.0)
        return 1.0;
    else
    {
        const double Re_t = std::pow(turbulent_kinetic_energy, 2) /
                            (kinematic_viscosity * turbulent_energy_dissipation_rate);
        const double f2 = 1.0 - 0.22 * std::exp(-1.0 * std::pow(Re_t / 6.0, 2));

        return f2;
    }
}

void CalculateStabilizationTau(double& tau,
                               double& element_length,
                               const array_1d<double, 3>& rVelocity,
                               const Matrix& rContravariantMetricTensor,
                               const double reaction,
                               const double effective_kinematic_viscosity,
                               const double alpha,
                               const double gamma,
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
    const double stab_reaction = std::pow(reaction, 2);

    PrintIfVariableIsNegative(stab_convection);

    const double velocity_norm = norm_2(rVelocity);
    tau = 1.0 / std::sqrt(stab_dynamics + stab_convection + stab_diffusion + stab_reaction);
    element_length = 2.0 * velocity_norm / std::sqrt(stab_convection);
}

double CalculateStabilizationTau(const double velocity_magnitude,
                                 const double length,
                                 const double effective_kinematic_viscosity,
                                 const double reaction,
                                 const double delta_time)
{
    const double stab_1 = 4.0 * effective_kinematic_viscosity / std::pow(length, 2);
    const double stab_2 = 2.0 * velocity_magnitude / length;
    const double stab_3 = reaction;
    const double stab_4 = 1.0 / delta_time;

    return 1.0 / (stab_1 + stab_2 + stab_3 + stab_4);
}

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy)
{
    // CheckIfVariableIsPositive(turbulent_kinematic_viscosity);
    // CheckIfVariableIsPositive(turbulent_kinetic_energy);

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

double CalculateGamma(const double C_mu,
                      const double f_mu,
                      const double turbulent_kinetic_energy,
                      const double turbulent_kinematic_viscosity)
{
    return std::max<double>(
        0.0, C_mu * f_mu * turbulent_kinetic_energy / turbulent_kinematic_viscosity);
}

void CalculateCrossWindDiffusionParameters(double& chi,
                                           double& k1,
                                           double& k2,
                                           const double velocity_magnitude,
                                           const double tau,
                                           const double effective_kinematic_viscosity,
                                           const double reaction,
                                           const double alpha,
                                           const double gamma,
                                           const double delta_time,
                                           const double element_length)
{
    const double reaction_dynamics = reaction + (1 - alpha) / (gamma * delta_time);

    chi = 2.0 / (std::abs(reaction_dynamics) * element_length + 2.0 * velocity_magnitude);

    double value = 0.0;
    value =
        std::abs(0.5 * (velocity_magnitude - tau * velocity_magnitude * reaction_dynamics +
                        tau * velocity_magnitude * std::abs(reaction_dynamics))) *
            element_length -
        (effective_kinematic_viscosity + tau * std::pow(velocity_magnitude, 2)) +
        (reaction_dynamics + tau * reaction_dynamics * std::abs(reaction_dynamics)) *
            std::pow(element_length, 2) / 6.0;

    k1 = std::max(value, 0.0);

    value = std::abs(0.5 * (velocity_magnitude + tau * velocity_magnitude *
                                                     std::abs(reaction_dynamics))) *
                element_length -
            effective_kinematic_viscosity +
            (reaction_dynamics + tau * reaction_dynamics * std::abs(reaction_dynamics)) *
                std::pow(element_length, 2) / 6.0;
    k2 = std::max(value, 0.0);
}

void CalculatePositiveValuesList(Vector& rOutput, const Vector& rInput)
{
    const int n = rInput.size();
    rOutput.resize(n);

    for (int i = 0; i < n; ++i)
    {
        if (rInput[i] >= std::numeric_limits<double>::epsilon())
            rOutput[i] = 1.0;
        else
            rOutput[i] = 0.0;
    }
}

void CalculateTurbulentValues(double& turbulent_kinetic_energy,
                              double& turbulent_energy_dissipation_rate,
                              const double y_plus,
                              const double kinematic_viscosity,
                              const double wall_distance,
                              const double c_mu,
                              const double von_karman)
{
    const double u_tau = y_plus * kinematic_viscosity / wall_distance;
    turbulent_kinetic_energy = std::pow(u_tau, 2) / std::sqrt(c_mu);
    turbulent_energy_dissipation_rate =
        std::pow(u_tau, 3) / (von_karman * wall_distance);
}

void CalculateTurbulentValues(double& turbulent_kinetic_energy,
                              double& turbulent_energy_dissipation_rate,
                              const double velocity_mag,
                              const double turbulence_intensity,
                              const double mixing_length,
                              const double c_mu)
{
    turbulent_kinetic_energy = 1.5 * std::pow(velocity_mag * turbulence_intensity, 2);
    turbulent_energy_dissipation_rate = c_mu * std::pow(turbulent_kinetic_energy, 1.5) / mixing_length;
}

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double, const double);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos
