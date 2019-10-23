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
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_UTILITIES)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_UTILITIES

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

namespace StabilizedConvectionDiffusionReactionAdjointUtilities
{
/**
 * @brief Calculates stabilization tau scalar derivatives
 *
 * \[
 *  -\tau_\phi^3\left[144\frac{\nu_\phi}{h^4_2}\left(\nu_{\phi,w}\right)^c + s_\phi\left(s_{\phi,w}\right)^c\right]
 * \]
 *
 * Where $w$ is the derivative variable
 *
 * @param rOutput                                        Scalar derivatives for each node w.r.t. $w$
 * @param Tau                                            Stabilization tau
 * @param EffectiveKinematicViscosity                    Effective kinematic viscosity $\nu_\phi$
 * @param Reaction                                       Reaction coefficient $s_\phi$
 * @param ElementLength                                  Element length $h_2$
 * @param rEffectiveKinematicViscosityScalarDerivatives  Scalar derivatives of effective kinematic viscosity $\left(\nu_{\phi,w}\right)$
 * @param rReactionScalarDerivatives                     Reaction scalar derivatives $\left(s_{\phi,w}\right)$
 */
inline void CalculateStabilizationTauScalarDerivatives(Vector& rOutput,
                                                       const double Tau,
                                                       const double EffectiveKinematicViscosity,
                                                       const double Reaction,
                                                       const double ElementLength,
                                                       const Vector& rEffectiveKinematicViscosityScalarDerivatives,
                                                       const Vector& rReactionScalarDerivatives)
{
    noalias(rOutput) =
        (rEffectiveKinematicViscosityScalarDerivatives *
             (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4)) +
         rReactionScalarDerivatives * (Reaction)) *
        (-1.0 * std::pow(Tau, 3));
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateCrossWindDiffusionCoeffVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double psi_one,
    const double element_length,
    const Matrix& rPsiOneDerivatives,
    const Matrix& rPsiTwoDerivatives,
    const Matrix& rEffectiveKinematicViscosityDerivatives,
    const Matrix& rElementLengthDerivatives)
{
    const double abs_psi_one = std::abs(psi_one);

    noalias(rOutput) =
        rPsiOneDerivatives * (0.5 * psi_one * element_length /
                              (abs_psi_one + std::numeric_limits<double>::epsilon())) +
        rElementLengthDerivatives * (0.5 * abs_psi_one) -
        rEffectiveKinematicViscosityDerivatives + rPsiTwoDerivatives;
}

inline double CalculateCrossWindDiffusionCoeffShapeSensitivity(
    const double psi_one,
    const double psi_one_deriv,
    const double element_length,
    const double element_length_deriv,
    const double effective_kinematic_viscosity_deriv,
    const double psi_two_deriv)
{
    const double abs_psi_one = std::abs(psi_one);
    return 0.5 * psi_one * element_length * psi_one_deriv /
               (abs_psi_one + std::numeric_limits<double>::epsilon()) +
           0.5 * abs_psi_one * element_length_deriv -
           effective_kinematic_viscosity_deriv + psi_two_deriv;
}

template <std::size_t TNumNodes>
inline void CalculateCrossWindDiffusionCoeffScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double psi_one,
    const double element_length,
    const Vector& rPsiOneScalarDerivatives,
    const Vector& rPsiTwoScalarDerivatives,
    const Vector& rEffectiveKinematicViscosityScalarDerivatives)
{
    noalias(rOutput) = rPsiOneScalarDerivatives *
                       (0.5 * psi_one * element_length / (std::abs(psi_one)) +
                        std::numeric_limits<double>::epsilon());
    noalias(rOutput) -= rEffectiveKinematicViscosityScalarDerivatives;
    noalias(rOutput) += rPsiTwoScalarDerivatives;
}

} // namespace StabilizedConvectionDiffusionReactionAdjointUtilities

} // namespace Kratos
#endif
