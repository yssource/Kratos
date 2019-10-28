//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//
// System includes
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"

#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_element.h"

#include "test_evm_k_epsilon_stabilized_cdr_adjoint_interface_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityVelocityDerivatives,
                          RansStabilizedCDRAdjointInterfaces)
{
    const unsigned int TNumNodes = 3;
    const unsigned int TDim = 2;
    using adjoint_element = RansEvmKAdjoint<TDim, TNumNodes>;
    using primal_element = RansEvmKElement<TDim, TNumNodes>;

    auto calculate_scalar_sensitivity =
        [TNumNodes, TDim](
            BoundedMatrix<double, TNumNodes, TDim>& rOutput,
            const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateEffectiveKinematicViscosityVelocityDerivatives(
                rOutput, rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunVectorKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", VELOCITY, calculate_scalar_sensitivity,
        calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateEffectiveKinematicViscosityShapeSensitivity,
                          RansStabilizedCDRAdjointInterfaces)
{
    const unsigned int TNumNodes = 3;
    const unsigned int TDim = 2;
    using adjoint_element = RansEvmKAdjoint<TDim, TNumNodes>;
    using primal_element = RansEvmKElement<TDim, TNumNodes>;

    auto calculate_scalar_sensitivity =
        [TNumNodes, TDim](
            BoundedMatrix<double, TNumNodes, TDim>& rOutput,
            const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateEffectiveKinematicViscosityVelocityDerivatives(
                rOutput, rData, rShapeFunctions, rShapeDerivatives, rProcessInfo);
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunVectorKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", VELOCITY, calculate_scalar_sensitivity,
        calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmKAdjoint<2, 3>;
    using primal_element = RansEvmKElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmK", "RansEvmKAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateEffectiveKinematicViscosityScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateEffectiveKinematicViscosityScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateEffectiveKinematicViscosity(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateReactionTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateReactionTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateReactionTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-6, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_KINETIC_ENERGY,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_KINETIC_ENERGY;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);
            return rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            return rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonAdjoint2D3N_CalculateSourceTermScalarDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          RansStabilizedCDRAdjointInterfaces)
{
    using adjoint_element = RansEvmEpsilonAdjoint<2, 3>;
    using primal_element = RansEvmEpsilonElement<2, 3>;
    const Variable<double>& perturb_variable = TURBULENT_ENERGY_DISSIPATION_RATE;

    auto calculate_scalar_sensitivity =
        [perturb_variable](
            BoundedVector<double, 3>& rOutput, const Vector& rShapeFunctions,
            const Matrix& rShapeDerivatives,
            const typename adjoint_element::BaseType& rElement,
            const typename adjoint_element::BaseType::ConvectionDiffusionReactionAdjointDataType& rData,
            const ProcessInfo& rProcessInfo, const int GaussIndex) {
            rElement.CalculateSourceTermScalarDerivatives(
                rOutput, perturb_variable, rData, rShapeFunctions,
                rShapeDerivatives, rProcessInfo);

            double temp = rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
            return temp;
        };

    auto calculate_scalar_value =
        [](const Vector& rShapeFunctions, const Matrix& rShapeDerivatives,
           const typename primal_element::BaseType& rElement,
           const typename primal_element::BaseType::ConvectionDiffusionReactionDataType& rData,
           const ProcessInfo& rProcessInfo, const int GaussIndex) {
            double temp = rElement.CalculateSourceTerm(
                rData, rShapeFunctions, rShapeDerivatives, rProcessInfo, 0);
            return temp;
        };

    RunScalarKEpsilon2D3NElementTest<primal_element, adjoint_element>(
        "RansEvmEpsilon", "RansEvmEpsilonAdjoint", perturb_variable,
        calculate_scalar_sensitivity, calculate_scalar_value, 1e-8, 1e-6, 1e-6);
}
} // namespace Testing
} // namespace Kratos