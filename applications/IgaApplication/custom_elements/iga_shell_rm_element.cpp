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
#include "iga_shell_rm_element.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaShellRMElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaShellRMElement>(NewId, geometry,
        pProperties);
}

void IgaShellRMElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
        SetDof(rElementalDofList, i, 3, SHEAR_A);
        SetDof(rElementalDofList, i, 4, SHEAR_B);
    }

    KRATOS_CATCH("")
}

void IgaShellRMElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetEquationId(rResult, i, 2, DISPLACEMENT_Z);
        SetEquationId(rResult, i, 3, SHEAR_A);
        SetEquationId(rResult, i, 4, SHEAR_B);
    }

    KRATOS_CATCH("")
}

void IgaShellRMElement::Initialize()
{
}

void IgaShellRMElement::CalculateAll(
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

void IgaShellRMElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaShellRMElement\" #" << Id();
}

} // namespace Kratos