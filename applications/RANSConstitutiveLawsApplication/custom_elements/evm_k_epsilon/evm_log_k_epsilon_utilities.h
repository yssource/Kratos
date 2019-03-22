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

#if !defined(KRATOS_EVM_LOG_K_EPSILON_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_LOG_K_EPSILON_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "custom_utilities/calculation_utilities.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "input_output/logger.h"

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
                                   const double f_mu);

double CalculateFmu(const double y_plus);

double CalculateF2(const double log_turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double log_turbulent_energy_dissipation_rate);

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity);

double CalculateGamma(const double log_turbulent_kinetic_energy,
                      const double log_turbulent_energy_dissipation_rate);

void CalculateStabilizationTau(double& tau,
                               double& element_length,
                               const array_1d<double, 3>& rVelocity,
                               const Matrix& rContravariantMetricTensor,
                               const double effective_kinematic_viscosity,
                               const double delta_time);

void CalculateCrossWindDiffusionParameters(double& chi,
                                           double& k1,
                                           double& k2,
                                           const double velocity_magnitude,
                                           const double tau,
                                           const double effective_kinematic_viscosity,
                                           const double reaction,
                                           const double element_length);

} // namespace EvmLogKepsilonModelUtilities

///@}

} // namespace Kratos

#endif