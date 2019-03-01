//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

#if !defined(KRATOS_EVM_K_EPSILON_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_K_EPSILON_UTILITIES_H_INCLUDED

#include "includes/define.h"
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
                                   const double f_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double mixing_length,
                                   const double minimum_viscosity);

double CalculateFmu(const double y_plus);

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate);

double CalculateFrictionVelocity(const double kinematic_viscosity,
                                 const double tangential_velocity_wall_gradient);

double CalculateYplus(const double friction_velocity,
                      const double wall_distance,
                      const double kinematic_viscosity);

double CalculateUTau(const double velocity_magnitude,
                     const double wall_distance,
                     const double kinematic_viscosity,
                     const double beta,
                     const double von_karman);

double CalculateStabilizationTau(const double velocity_magnitude,
                                 const double length,
                                 const double turbulent_kinetic_energy);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos

#endif