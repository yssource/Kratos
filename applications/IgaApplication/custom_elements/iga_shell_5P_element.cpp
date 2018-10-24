/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes
#include "includes/define.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "iga_shell_5P_element.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaShell5PElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaShell5PElement>(NewId, geometry,
        pProperties);
}

void IgaShell5PElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
        SetElementDof(rElementalDofList, i, 3, SHEAR_A);
        SetElementDof(rElementalDofList, i, 4, SHEAR_B);
    }

    KRATOS_CATCH("")
}

void IgaShell5PElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetElementEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetElementEquationId(rResult, i, 2, DISPLACEMENT_Z);
        SetElementEquationId(rResult, i, 3, SHEAR_A);
        SetElementEquationId(rResult, i, 4, SHEAR_B);
    }

    KRATOS_CATCH("")
}

void IgaShell5PElement::Initialize()
{
}

void IgaShell5PElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

    // get integration data

    const double integration_weight = GetValue(INTEGRATION_WEIGHT);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    const double E = GetValue(YOUNG_MODULUS);
    const double thickness = GetValue(THICKNESS);
    const double poisson_ratio = GetValue(POISSON_RATIO);

    KRATOS_CATCH("")
}

void IgaShell5PElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaShell5PElement\" #" << Id();
}

} // namespace Kratos