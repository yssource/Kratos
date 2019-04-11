//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_EVM_K_EPSILON_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_K_EPSILON_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"

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

///@}
///@name Kratos Classes
///@{

/**
 * @class VariableUtils
 * @ingroup KratosCore
 * @brief This class implements a set of auxiliar, already parallelized, methods to
 * perform some common tasks related with the variable values and fixity.
 * @details The methods are exported to python in order to add this improvements to the python interface
 * @author Riccardo Rossi
 * @author Ruben Zorrilla
 * @author Vicente Mataix Ferrandiz
 */

class EvmKepsilonModelUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// We create the Pointer related to EvmKepsilonModelUtilities
    KRATOS_CLASS_POINTER_DEFINITION(EvmKepsilonModelUtilities);

    template <unsigned int TDim>
    double static CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                                      const double turbulent_kinematic_viscosity,
                                      const double turbulent_kinetic_energy);

    double static CalculateTurbulentViscosity(const double C_mu,
                                              const double turbulent_kinetic_energy,
                                              const double turbulent_energy_dissipation_rate,
                                              const double f_mu);

    double static CalculateFmu(const double y_plus);

    double static CalculateF2(const double turbulent_kinetic_energy,
                              const double kinematic_viscosity,
                              const double turbulent_energy_dissipation_rate);

    double static CalculateGamma(const double C_mu,
                                 const double f_mu,
                                 const double turbulent_kinetic_energy,
                                 const double turbulent_kinematic_viscosity);

    void static CalculateTurbulentValues(double& turbulent_kinetic_energy,
                                         double& turbulent_energy_dissipation_rate,
                                         const double y_plus,
                                         const double kinematic_viscosity,
                                         const double wall_distance,
                                         const double c_mu,
                                         const double von_karman);

    void static CalculateTurbulentValues(double& turbulent_kinetic_energy,
                                         double& turbulent_energy_dissipation_rate,
                                         const double velocity_mag,
                                         const double turbulence_intensity,
                                         const double mixing_length,
                                         const double c_mu);

    void static CalculatePositiveValuesList(Vector& rOutput, const Vector& rInput);

    void CalculateTurbulentViscosityForModelPart(ModelPart& rModelPart,
                                                 unsigned int Step = 0);

    void UpdateBoundaryConditions(ModelPart& rModelPart);

    void AssignInitialValues(ModelPart& rModelPart);

}; // class EvmKepsilonModelUtilities

///@}

} // namespace Kratos

#endif