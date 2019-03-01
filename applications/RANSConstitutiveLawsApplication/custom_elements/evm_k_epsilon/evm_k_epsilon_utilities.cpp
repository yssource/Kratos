
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
                                   const double f_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double mixing_length,
                                   const double minimum_viscosity)
{
    KRATOS_DEBUG_ERROR_IF(turbulent_kinetic_energy < 0.0)
        << "TURBULENT_KINETIC_ENERGY < 0.0 [ " << std::scientific
        << turbulent_kinetic_energy << " < 0.0 ] in CalculateTurbulentViscosity\n";

    KRATOS_DEBUG_ERROR_IF(turbulent_energy_dissipation_rate < 0.0)
        << "TURBULENT_ENERGY_DISSIPATION_RATE < 0.0 [ " << std::scientific
        << turbulent_energy_dissipation_rate
        << " < 0.0 ] in CalculateTurbulentViscosity\n";

    const double limited_mixing_length = std::min<double>(
        C_mu * std::pow(turbulent_kinetic_energy, 1.5) / turbulent_energy_dissipation_rate,
        mixing_length);

    const double nu_t = std::max<double>(
        minimum_viscosity, f_mu * limited_mixing_length * std::sqrt(turbulent_kinetic_energy));

    return nu_t;
}

double CalculateFmu(const double y_plus)
{
    KRATOS_DEBUG_ERROR_IF(y_plus < 0.0) << "Y_PLUS < 0.0 [ " << std::scientific
                                        << y_plus << " < 0.0 ] in CalculateFmu\n";

    return 1.0 - std::exp(-0.0115 * y_plus);
}

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate)
{
    KRATOS_DEBUG_ERROR_IF(turbulent_kinetic_energy < 0.0)
        << "TURBULENT_KINETIC_ENERGY < 0.0 [ " << std::scientific
        << turbulent_kinetic_energy << " < 0.0 ] in CalculateF2\n";

    KRATOS_DEBUG_ERROR_IF(turbulent_energy_dissipation_rate < 0.0)
        << "TURBULENT_ENERGY_DISSIPATION_RATE < 0.0 [ " << std::scientific
        << turbulent_energy_dissipation_rate << " < 0.0 ] in CalculateF2\n";

    KRATOS_DEBUG_ERROR_IF(kinematic_viscosity < 0.0)
        << "KINEMATIC_VISCOSITY < 0.0 [ " << std::scientific
        << kinematic_viscosity << " < 0.0 ] in CalculateF2\n";

    return 1.0 - 0.22 * std::exp(-1.0 * std::pow(std::pow(turbulent_kinetic_energy, 2) /
                                                     (6 * kinematic_viscosity * turbulent_energy_dissipation_rate),
                                                 2));
}

double CalculateFrictionVelocity(const double kinematic_viscosity,
                                 const double tangential_velocity_wall_gradient)
{
    KRATOS_DEBUG_ERROR_IF(kinematic_viscosity < 0.0)
        << "KINEMATIC_VISCOSITY < 0.0 [ " << std::scientific
        << kinematic_viscosity << " < 0.0 ] in CalculateFrictionVelocity\n";

    KRATOS_DEBUG_ERROR_IF(tangential_velocity_wall_gradient < 0.0)
        << "Tangential_wall_velocity_gradient < 0.0 [ " << std::scientific
        << tangential_velocity_wall_gradient << " < 0.0 ] in CalculateFrictionVelocity\n";

    return std::sqrt(kinematic_viscosity * tangential_velocity_wall_gradient);
}

double CalculateYplus(const double friction_velocity, const double wall_distance, const double kinematic_viscosity)
{
    KRATOS_DEBUG_ERROR_IF(kinematic_viscosity < 0.0)
        << "KINEMATIC_VISCOSITY < 0.0 [ " << std::scientific
        << kinematic_viscosity << " < 0.0 ] in CalculateYplus\n";

    KRATOS_DEBUG_ERROR_IF(friction_velocity < 0.0)
        << "friction_velocity < 0.0 [ " << std::scientific << friction_velocity
        << " < 0.0 ] in CalculateYplus\n";

    KRATOS_DEBUG_ERROR_IF(wall_distance < 0.0)
        << "wall_distance < 0.0 [ " << std::scientific << wall_distance
        << " < 0.0 ] in CalculateYplus\n";

    return friction_velocity * wall_distance / kinematic_viscosity;
}

double CalculateUTau(const double velocity_magnitude,
                     const double wall_distance,
                     const double kinematic_viscosity,
                     const double beta,
                     const double von_karman)
{
    KRATOS_DEBUG_ERROR_IF(kinematic_viscosity < 0.0)
        << "KINEMATIC_VISCOSITY < 0.0 [ " << std::scientific
        << kinematic_viscosity << " < 0.0 ] in CalculateUTau\n";

    KRATOS_DEBUG_ERROR_IF(wall_distance < 0.0)
        << "wall_distance < 0.0 [ " << std::scientific << wall_distance
        << " < 0.0 ] in CalculateUTau\n";

    KRATOS_DEBUG_ERROR_IF(velocity_magnitude < 0.0)
        << "velocity_magnitude < 0.0 [ " << std::scientific
        << velocity_magnitude << " < 0.0 ] in CalculateUTau\n";

    const unsigned int max_iterations = 10;
    double u_tau = kinematic_viscosity / wall_distance;
    for (unsigned int i = 0; i < max_iterations; ++i)
    {
        u_tau = kinematic_viscosity *
                std::exp(((velocity_magnitude / u_tau) - beta) * von_karman) / wall_distance;
    }

    double y_plus = u_tau * wall_distance / kinematic_viscosity;

    if (y_plus < 11.06)
        u_tau = std::sqrt(velocity_magnitude * kinematic_viscosity / wall_distance);

    KRATOS_ERROR_IF(u_tau < 0.0) << "Calculated u_tau < 0.0 [ "
                                 << std::scientific << u_tau << " < 0.0 ]\n";

    return u_tau;
}

double CalculateStabilizationTau(const double velocity_magnitude,
                                 const double length,
                                 const double turbulent_kinetic_energy)
{
    double alpha(1.0), Pe(0.0);
    if (turbulent_kinetic_energy != 0.0)
        Pe = std::max<double>(velocity_magnitude * length / (2.0 * turbulent_kinetic_energy), 0.0);
    if (Pe != 0.0)
        alpha = (std::exp(2.0 * Pe) + 1) / (std::exp(2.0 * Pe) - 1) - 1.0 / Pe;

    if (velocity_magnitude != 0.0)
        return alpha * length / (2.0 * velocity_magnitude);
    else
        return 0.0;
}

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos
