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
#include "includes/ublas_interface.h"
#include "input_output/logger.h"
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
#define CheckIfVariableIsPositive(variable) \
    KRATOS_DEBUG_ERROR_IF(variable < 0.0)   \
        << #variable << " < 0.0 [ " << std::scientific << variable << " < 0.0 ]\n";

double CalculateTurbulentViscosity(const double C_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double minimum_viscosity);

double CalculateFmu(const double y_plus);

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate);

double CalculateYplus(const double velocity_norm,
                      const double wall_distance,
                      const double kinematic_viscosity,
                      const double von_karman,
                      const double beta,
                      const unsigned int max_iterations);

double CalculateStabilizationTau(const double velocity_magnitude,
                                 const double length,
                                 const double effective_kinematic_viscosity,
                                 const double delta_time);

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos

#endif