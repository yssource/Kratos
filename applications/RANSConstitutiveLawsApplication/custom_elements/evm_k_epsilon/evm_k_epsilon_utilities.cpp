
#include "evm_k_epsilon_utilities.h"

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
                                   const double minimum_viscosity)
{
    // CheckIfVariableIsPositive(turbulent_kinetic_energy);
    // CheckIfVariableIsPositive(turbulent_energy_dissipation_rate);

    // const double limited_mixing_length = std::min<double>(
    //     C_mu * std::pow(turbulent_kinetic_energy, 1.5) /
    //     turbulent_energy_dissipation_rate, mixing_length);

    // const double nu_t = std::max<double>(
    //     minimum_viscosity, f_mu * limited_mixing_length * std::sqrt(turbulent_kinetic_energy));

    return std::max<double>(minimum_viscosity,
                            C_mu * std::pow(turbulent_kinetic_energy, 2) /
                                turbulent_energy_dissipation_rate);
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

    return 1.0 - 0.22 * std::exp(-1.0 * std::pow(std::pow(turbulent_kinetic_energy, 2) /
                                                     (6 * kinematic_viscosity * turbulent_energy_dissipation_rate),
                                                 2));
}

double CalculateYplus(const double velocity_norm,
                      const double wall_distance,
                      const double kinematic_viscosity,
                      const double von_karman,
                      const double beta,
                      const unsigned int max_iterations)
{
    CheckIfVariableIsPositive(velocity_norm);
    CheckIfVariableIsPositive(wall_distance);
    CheckIfVariableIsPositive(kinematic_viscosity);
    CheckIfVariableIsPositive(von_karman);
    CheckIfVariableIsPositive(beta);

    // try linear law
    double y_plus = std::sqrt(velocity_norm * wall_distance / kinematic_viscosity);

    // If the linear low doesnt match within the range, try logrithmic law
    if (y_plus > 11.06)
    {
        unsigned int i;
        double u_tau = std::sqrt(velocity_norm * kinematic_viscosity / wall_distance);
        double prev_u_tau = 0.0;
        for (i = 0; i < max_iterations; ++i)
        {
            prev_u_tau = u_tau;
            u_tau = velocity_norm /
                    (std::log(u_tau * wall_distance / kinematic_viscosity) / von_karman + beta);
        }
        const double delta_u_tau = std::abs(u_tau - prev_u_tau);
        KRATOS_INFO_IF("TurbulenceEvmProcess", delta_u_tau > 1e-5)
            << "WARNING: Maximum number of iterations reached for y_plus "
               "calculation. error_u_tau = "
            << std::scientific << delta_u_tau << ".\n";
        y_plus = u_tau * wall_distance / kinematic_viscosity;
    }
    return y_plus;
}

double CalculateStabilizationTau(const double velocity_magnitude,
                                 const double length,
                                 const double effective_kinematic_viscosity,
                                 const double delta_time)
{
    const double Pe = 0.5 * velocity_magnitude * length / effective_kinematic_viscosity;
    const double alpha = Pe * (std::exp(2.0 * Pe) + 1) / (std::exp(2.0 * Pe) - 1) - 1.0;

    const double stab_1 =
        std::pow(velocity_magnitude, 2) / (effective_kinematic_viscosity * alpha);
    const double stab_2 = effective_kinematic_viscosity / std::pow(length, 2);
    const double stab_3 = 1.0 / delta_time;

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

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos
