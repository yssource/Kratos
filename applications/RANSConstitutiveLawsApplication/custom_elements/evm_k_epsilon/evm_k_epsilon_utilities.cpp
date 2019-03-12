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
                                   const double f_mu,
                                   const double minimum_viscosity,
                                   const double maximum_viscosity)
{
    CheckIfVariableIsPositive(turbulent_kinetic_energy);
    CheckIfVariableIsPositive(turbulent_energy_dissipation_rate);

    const double nu_t = std::min<double>(
        std::max<double>(minimum_viscosity, C_mu * f_mu * std::pow(turbulent_kinetic_energy, 2) /
                                                turbulent_energy_dissipation_rate),
        maximum_viscosity);

    return nu_t;
}

double CalculateFmu(const double y_plus)
{
    CheckIfVariableIsPositive(y_plus);
    return 1.0 - std::exp(-0.0115 * y_plus);
}

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate)
{
    CheckIfVariableIsPositive(turbulent_kinetic_energy);
    CheckIfVariableIsPositive(turbulent_energy_dissipation_rate);
    CheckIfVariableIsPositive(kinematic_viscosity);

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

double CalculateStabilizationTau(const double velocity_magnitude,
                                 const double length,
                                 const double effective_kinematic_viscosity,
                                 const double reaction,
                                 const double delta_time)
{
    CheckIfVariableIsPositive(velocity_magnitude);
    CheckIfVariableIsPositive(length);
    CheckIfVariableIsPositive(effective_kinematic_viscosity);
    CheckIfVariableIsPositive(delta_time);

    PrintIfVariableIsNegative(velocity_magnitude);
    PrintIfVariableIsNegative(length);
    PrintIfVariableIsNegative(effective_kinematic_viscosity);
    PrintIfVariableIsNegative(reaction);
    PrintIfVariableIsNegative(delta_time);

    // const double Pe = std::max<double>(
    //     0.5 * velocity_magnitude * length / effective_kinematic_viscosity, 1e-15);
    // const double alpha = std::max<double>(
    //     std::min<double>((std::exp(2.0 * Pe) + 1) / (std::exp(2.0 * Pe) - 1) - 1.0 / Pe, 1.0),
    //     0.0);

    // const double stab_tau = alpha * length / (2 * velocity_magnitude);

    // PrintIfVariableIsNegative(stab_tau);
    // return stab_tau;

    const double stab_1 = 4.0 * effective_kinematic_viscosity / std::pow(length, 2);
    const double stab_2 = 2.0 * velocity_magnitude / length;
    const double stab_3 = reaction;

    return 1.0 / (stab_1 + stab_2 + stab_3);
}

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient)
{
    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            source += symmetric_velocity_gradient(i, j) *
                      symmetric_velocity_gradient(i, j);

    return 0.5 * source;
}

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy)
{
    CheckIfVariableIsPositive(turbulent_kinematic_viscosity);
    CheckIfVariableIsPositive(turbulent_kinetic_energy);

    const double velocity_divergence =
        CalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) =
        turbulent_kinematic_viscosity *
            (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity) -
        (2.0 / 3.0) * turbulent_kinetic_energy * identity;

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

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&);
template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double, const double);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos
