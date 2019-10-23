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

inline double CalculateScalarProduct(const Vector& rVector1, const array_1d<double, 3>& rVector2)
{
    double result = 0.0;
    for (std::size_t i_dim = 0; i_dim < rVector1.size(); ++i_dim)
        result += rVector1[i_dim] * rVector2[i_dim];
    return result;
}

inline void CalculateAbsoluteScalarValueScalarDerivatives(Vector& rOutput,
                                                          const double scalar_value,
                                                          const Vector& rScalarValueDerivatives)
{
    noalias(rOutput) =
        rScalarValueDerivatives *
        (scalar_value / (std::abs(scalar_value) + std::numeric_limits<double>::epsilon()));
}

template <std::size_t TNumNodes>
inline void CalculateAbsoluteScalarValueScalarDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                                          const double scalar_value,
                                                          const Vector& rScalarValueDerivatives)
{
    noalias(rOutput) =
        rScalarValueDerivatives *
        (scalar_value / (std::abs(scalar_value) + std::numeric_limits<double>::epsilon()));
}

template <std::size_t TNumNodes>
inline void CalculateAbsoluteScalarValueScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double scalar_value,
    const BoundedVector<double, TNumNodes>& rScalarValueDerivatives)
{
    noalias(rOutput) =
        rScalarValueDerivatives *
        (scalar_value / (std::abs(scalar_value) + std::numeric_limits<double>::epsilon()));
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateAbsoluteScalarValueVectorDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double scalar_value,
    const BoundedMatrix<double, TNumNodes, TDim>& rScalarValueDerivatives)
{
    noalias(rOutput) =
        rScalarValueDerivatives *
        (scalar_value / (std::abs(scalar_value) + std::numeric_limits<double>::epsilon()));
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateAbsoluteScalarValueVectorDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double scalar_value,
    const Matrix& rScalarValueDerivatives)
{
    noalias(rOutput) =
        rScalarValueDerivatives *
        (scalar_value / (std::abs(scalar_value) + std::numeric_limits<double>::epsilon()));
}

inline void CalculateAbsoluteScalarValueVectorDerivatives(Matrix& rOutput,
                                                          const double scalar_value,
                                                          const Matrix& rScalarValueDerivatives)
{
    noalias(rOutput) =
        rScalarValueDerivatives *
        (scalar_value / (std::abs(scalar_value) + std::numeric_limits<double>::epsilon()));
}

template<std::size_t TNumNodes>
inline void CalculatePsiOneScalarDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                      const double velocity_norm,
                                      const double reaction_tilde,
                                      const double tau,
                                      const Vector& rTauScalarDerivatives,
                                      const Vector& rAbsoluteReactionTildeScalarDerivatives)
{
    const double absolute_reaction_tilde = std::abs(reaction_tilde);

    noalias(rOutput) = rTauScalarDerivatives * (velocity_norm * absolute_reaction_tilde);
    noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * velocity_norm);
}

template<std::size_t TDim, std::size_t TNumNodes>
inline void CalculatePsiOneVelocityDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                        const double velocity_norm,
                                        const double reaction_tilde,
                                        const double tau,
                                        const Matrix& rTauDerivatives,
                                        const Matrix& rAbsoluteReactionTildeDerivatives,
                                        const Matrix& rVelocityMagnitudeDerivatives)
{
    noalias(rOutput) = rVelocityMagnitudeDerivatives +
                       rTauDerivatives * (velocity_norm * reaction_tilde) +
                       rVelocityMagnitudeDerivatives * (tau * reaction_tilde) +
                       rAbsoluteReactionTildeDerivatives * (tau * velocity_norm);
}

inline double CalculatePsiOneShapeSensitivity(const double tau,
                                       const double tau_deriv,
                                       const double velocity_magnitude,
                                       const double reaction,
                                       const double reaction_deriv,
                                       const double bossak_alpha,
                                       const double bossak_gamma,
                                       const double delta_time,
                                       const double DynamicTau)
{
    const double reaction_dynamics =
        reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
    const double abs_reaction_dynamics = std::abs(reaction_dynamics);

    return tau_deriv * velocity_magnitude * abs_reaction_dynamics +
           tau * velocity_magnitude * reaction_dynamics * reaction_deriv /
               (abs_reaction_dynamics + std::numeric_limits<double>::epsilon());
}

template <std::size_t TNumNodes>
inline void CalculateChiScalarDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                          const double Chi,
                                          const double ElementLength,
                                          const double BossakAlpha,
                                          const double BossakGamma,
                                          const double DeltaTime,
                                          const double Reaction,
                                          const double DynamicTau,
                                          const Vector& rReactionScalarDerivatives)
{
    const double reaction_tilde =
        Reaction + DynamicTau * (1 - BossakAlpha) / (BossakGamma * DeltaTime);

    CalculateAbsoluteScalarValueScalarDerivatives(rOutput, reaction_tilde,
                                                  rReactionScalarDerivatives);
    noalias(rOutput) = rOutput * (-0.5 * std::pow(Chi, 2) * ElementLength);
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateChiVelocityDerivatives(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                            const double Chi,
                                            const double ElementLength,
                                            const double BossakAlpha,
                                            const double BossakGamma,
                                            const double DeltaTime,
                                            const double Reaction,
                                            const double DynamicTau,
                                            const Matrix& rReactionDerivatives,
                                            const Matrix& rVelocityMagnitudeDerivatives,
                                            const Matrix& rElementLengthDerivatives)
{
    const double reaction_tilde =
        Reaction + DynamicTau * (1 - BossakAlpha) / (BossakGamma * DeltaTime);
    const double abs_reaction_tilde = std::abs(reaction_tilde);

    CalculateAbsoluteScalarValueVectorDerivatives(rOutput, reaction_tilde, rReactionDerivatives);

    noalias(rOutput) = (rOutput * ElementLength + rElementLengthDerivatives * abs_reaction_tilde +
                        rVelocityMagnitudeDerivatives * 2.0) *
                       (-0.5 * std::pow(Chi, 2));
}

inline double CalculateChiShapeSensitivity(const double chi,
                                           const double reaction,
                                           const double reaction_deriv,
                                           const double element_length,
                                           const double element_length_deriv,
                                           const double bossak_alpha,
                                           const double bossak_gamma,
                                           const double delta_time,
                                           const double DynamicTau)
{
    const double reaction_tilde =
        reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
    const double abs_reaction_tilde = std::abs(reaction_tilde);

    return -0.5 * std::pow(chi, 2) *
           (abs_reaction_tilde * element_length_deriv +
            reaction_tilde * element_length * reaction_deriv /
                (abs_reaction_tilde + std::numeric_limits<double>::epsilon()));
}

template <std::size_t TNumNodes>
inline void CalculateResidualScalarDerivative(BoundedVector<double, TNumNodes>& rOutput,
                                              const double scalar_value,
                                              const double reaction,
                                              const array_1d<double, 3>& rVelocity,
                                              const Vector& rReactionScalarDerivatives,
                                              const Vector& rSourceScalarDerivatives,
                                              const Vector& rShapeFunctions,
                                              const Matrix& rShapeFunctionDerivatives,
                                              const Variable<double>& rPrimalVariable,
                                              const Variable<double>& rDerivativeVariable)
{
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
        double value = 0.0;

        value += scalar_value * rReactionScalarDerivatives[i_node];
        value -= rSourceScalarDerivatives[i_node];

        if (rPrimalVariable == rDerivativeVariable)
        {
            value += reaction * rShapeFunctions[i_node];
            value += CalculateScalarProduct(shape_function_gradient, rVelocity);
        }

        rOutput[i_node] = value;
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateResidualVelocityDerivative(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                                const double primal_variable_value,
                                                const array_1d<double, 3>& rPrimalVariableGradient,
                                                const Matrix& rReactionDerivatives,
                                                const Matrix& rSourceDerivatives,
                                                const Vector& rGaussShapeFunctions)
{
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
            rOutput(i_node, i_dim) =
                rGaussShapeFunctions[i_node] * rPrimalVariableGradient[i_dim];

    noalias(rOutput) = rOutput + rReactionDerivatives * primal_variable_value - rSourceDerivatives;
}

inline double CalculateResidualShapeSensitivity(const double residual,
                                                const array_1d<double, 3>& rVelocity,
                                                const Matrix& rShapeFunctionDerivShapeSensitivity,
                                                const double scalar_value,
                                                const Vector& rNodalScalarValues,
                                                const double reaction_deriv,
                                                const double source_deriv)
{
    const double abs_residual = std::abs(residual);

    if (abs_residual <= std::numeric_limits<double>::epsilon())
    {
        return 0.0;
    }
    else
    {
        Vector r_velocity(rShapeFunctionDerivShapeSensitivity.size2());
        for (std::size_t i = 0; i < rShapeFunctionDerivShapeSensitivity.size2(); ++i)
            r_velocity[i] = rVelocity[i];
        // const Vector& r_velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);
        Vector primal_variable_gradient_shape_sensitivity(
            rShapeFunctionDerivShapeSensitivity.size2());
        noalias(primal_variable_gradient_shape_sensitivity) =
            prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

        return residual *
               (inner_prod(r_velocity, primal_variable_gradient_shape_sensitivity) +
                reaction_deriv * scalar_value - source_deriv) /
               abs_residual;
    }
}

template <std::size_t TNumNodes>
inline void CalculatePositivityPreservationCoefficientScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double chi,
    const double residual,
    const double scalar_gradient_norm,
    const double velocity_norm_square,
    const Vector& rChiScalarDerivatives,
    const Vector& rAbsoluteResidualScalarDerivatives,
    const Vector& rAbsoluteScalarGradientScalarDerivative,
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rDerivativeVariable)
{
    const double abs_residual = std::abs(residual);

    if (scalar_gradient_norm <= std::numeric_limits<double>::epsilon() ||
        velocity_norm_square <= std::numeric_limits<double>::epsilon())
    {
        rOutput.clear();
    }
    else
    {
        noalias(rOutput) = rAbsoluteResidualScalarDerivatives *
                           (chi / (velocity_norm_square * scalar_gradient_norm));
        noalias(rOutput) +=
            rChiScalarDerivatives *
            (abs_residual / (velocity_norm_square * scalar_gradient_norm));

        if (rPrimalVariable == rDerivativeVariable)
            noalias(rOutput) -=
                rAbsoluteScalarGradientScalarDerivative *
                (chi * abs_residual / (std::pow(scalar_gradient_norm, 2) * velocity_norm_square));
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculatePositivityPreservationCoefficientVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double absolute_residual,
    const double primal_variable_gradient_norm,
    const double velocity_magnitude,
    const double chi,
    const Matrix& rChiDerivatives,
    const Matrix& rAbsoluteResidualDerivatives,
    const Matrix& rVelocityMagnitudeDerivatives)
{
    const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

    if (primal_variable_gradient_norm <= std::numeric_limits<double>::epsilon() ||
        velocity_magnitude_square <= std::numeric_limits<double>::epsilon())
    {
        rOutput.clear();
    }

    else
    {
        noalias(rOutput) =
            (rVelocityMagnitudeDerivatives * (-2.0 * chi / velocity_magnitude) + rChiDerivatives) *
                (absolute_residual / (velocity_magnitude_square * primal_variable_gradient_norm)) +
            rAbsoluteResidualDerivatives *
                (chi / (primal_variable_gradient_norm * velocity_magnitude_square));
    }
}

inline double CalculatePositivityPreservationCoefficientShapeSensitivity(
    const double chi,
    const double chi_deriv,
    const double abs_residual,
    const double abs_residual_deriv,
    const double velocity_magnitude_square,
    const double scalar_gradient_norm,
    const double scalar_gradient_norm_deriv)
{
    if (scalar_gradient_norm <= std::numeric_limits<double>::epsilon() ||
        velocity_magnitude_square <= std::numeric_limits<double>::epsilon())
    {
        return 0.0;
    }
    else
    {
        return chi_deriv * abs_residual / (scalar_gradient_norm * velocity_magnitude_square) +
               chi * abs_residual_deriv / (scalar_gradient_norm * velocity_magnitude_square) -
               chi * abs_residual * scalar_gradient_norm_deriv /
                   (std::pow(scalar_gradient_norm, 2) * velocity_magnitude_square);
    }
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

template <std::size_t TNumNodes>
inline void CalculateStreamLineDiffusionCoeffScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double element_length,
    const double tau,
    const double velocity_norm,
    const double reaction_tilde,
    const double psi_one,
    const double psi_two,
    const Vector& rPsiOneScalarDerivatives,
    const Vector& rPsiTwoScalarDerivatives,
    const Vector& rTauScalarDerivatives,
    const Vector& rReactionTildeScalarDerivatives,
    const Vector& rEffectiveViscosityScalarDerivatives)
{
    noalias(rOutput) = rPsiOneScalarDerivatives;
    noalias(rOutput) -= rTauScalarDerivatives * (velocity_norm * reaction_tilde);
    noalias(rOutput) -= rReactionTildeScalarDerivatives * (tau * velocity_norm);

    const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
    noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / (std::abs(coeff)) +
                                  std::numeric_limits<double>::epsilon());

    noalias(rOutput) += rPsiTwoScalarDerivatives;
    noalias(rOutput) -= rEffectiveViscosityScalarDerivatives;
    noalias(rOutput) -= rTauScalarDerivatives * std::pow(velocity_norm, 2);
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateStreamLineDiffusionCoeffVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double element_length,
    const double tau,
    const double velocity_norm,
    const double reaction_tilde,
    const double psi_one,
    const double psi_two,
    const Matrix& rVelocityMagnitudeDerivatives,
    const Matrix& rPsiOneDerivatives,
    const Matrix& rPsiTwoDerivatives,
    const Matrix& rTauDerivatives,
    const Matrix& rReactionTildeDerivatives,
    const Matrix& rEffectiveViscosityDerivatives,
    const Matrix& rElementLengthDerivatives)
{
    noalias(rOutput) =
        (rPsiOneDerivatives - rTauDerivatives * (velocity_norm * reaction_tilde) -
         rVelocityMagnitudeDerivatives * (tau * reaction_tilde) -
         rReactionTildeDerivatives * (tau * velocity_norm));

    const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
    noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / (std::abs(coeff)) +
                                  std::numeric_limits<double>::epsilon());

    noalias(rOutput) += rElementLengthDerivatives * (0.5 * std::abs(coeff));

    noalias(rOutput) += rPsiTwoDerivatives;
    noalias(rOutput) -= rEffectiveViscosityDerivatives;
    noalias(rOutput) -= rTauDerivatives * std::pow(velocity_norm, 2);
    noalias(rOutput) -= rVelocityMagnitudeDerivatives * (2.0 * tau * velocity_norm);
}

inline double CalculateStreamLineDiffusionCoeffShapeSensitivity(
    const double psi_one,
    const double psi_one_deriv,
    const double tau,
    const double tau_deriv,
    const double velocity_magnitude,
    const double reaction,
    const double reaction_deriv,
    const double element_length,
    const double element_length_deriv,
    const double effective_kinematic_viscosity_deriv,
    const double psi_two_deriv,
    const double bossak_alpha,
    const double bossak_gamma,
    const double delta_time,
    const double DynamicTau)
{
    const double reaction_dynamics =
        reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
    const double coeff = psi_one - tau * velocity_magnitude * reaction_dynamics;
    const double abs_coeff = std::abs(coeff);
    double shape_sensitivity = 0.0;

    shape_sensitivity += psi_one_deriv - tau_deriv * velocity_magnitude * reaction_dynamics -
                         tau * velocity_magnitude * reaction_deriv;
    shape_sensitivity *= 0.5 * coeff * element_length /
                         (abs_coeff + std::numeric_limits<double>::epsilon());
    shape_sensitivity += 0.5 * abs_coeff * element_length_deriv;
    shape_sensitivity -= effective_kinematic_viscosity_deriv +
                         tau_deriv * std::pow(velocity_magnitude, 2);
    shape_sensitivity += psi_two_deriv;

    return shape_sensitivity;
}

} // namespace StabilizedConvectionDiffusionReactionAdjointUtilities

} // namespace Kratos
#endif
