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
#include "iga_shell_3P_element.h"
#include "iga_application_variables.h"
#include "iga_debug.h"

namespace Kratos {

Element::Pointer IgaShell3PElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaShell3PElement>(NewId, geometry,
        pProperties);
}

void IgaShell3PElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void IgaShell3PElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetElementEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetElementEquationId(rResult, i, 2, DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void IgaShell3PElement::Initialize()
{
    KRATOS_TRY;

    KRATOS_CATCH("")
}

void IgaShell3PElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

    // temporary debug data
    auto expected_data = Parameters(GetValue(DEBUG_EXPECTED_DATA));

    // get integration data

    const double integration_weight = GetValue(INTEGRATION_WEIGHT);
    Vector& shape_function_values = GetValue(SHAPE_FUNCTION_VALUES);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    const auto& properties = GetProperties();

    const double E = properties.GetValue(YOUNG_MODULUS);
    const double thickness = properties.GetValue(THICKNESS);
    const double poisson_ratio = properties.GetValue(POISSON_RATIO);

    // material matrices

    BoundedMatrix<double, 3, 3> db;
    db.clear();
    db(0, 0) = 1.0;
    db(0, 1) = poisson_ratio;
    db(1, 0) = poisson_ratio;
    db(1, 1) = 1.0;
    db(2, 2) = (1.0 - poisson_ratio) / 2.0;
    db *= E * pow(thickness, 3) / (12.0 * (1.0 - poisson_ratio * poisson_ratio));

    BoundedMatrix<double, 3, 3> dm;
    dm.clear();
    dm(0, 0) = 1.0;
    dm(0, 1) = poisson_ratio;
    dm(1, 0) = poisson_ratio;
    dm(1, 1) = 1.0;
    dm(2, 2) = (1.0 - poisson_ratio) / 2.0;
    dm *= E * thickness / (1.0 - poisson_ratio * poisson_ratio);

    IgaDebug::CheckMatrix(expected_data, "_Db", db);
    IgaDebug::CheckMatrix(expected_data, "_Dm", dm);

    // reference and actual configuration

    Configuration ref;
    Configuration act;

    ref.a1.clear();
    ref.a2.clear();

    act.a1.clear();
    act.a2.clear();

    ref.a11.clear();
    ref.a12.clear();
    ref.a22.clear();

    act.a11.clear();
    act.a12.clear();
    act.a22.clear();

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        const Node<3>& node = GetGeometry()[i];
        ref.a1 += shape_derivatives(0, i) * node.GetInitialPosition();
        ref.a2 += shape_derivatives(1, i) * node.GetInitialPosition();
        ref.a11 += shape_derivatives(2, i) * node.GetInitialPosition();
        ref.a12 += shape_derivatives(3, i) * node.GetInitialPosition();
        ref.a22 += shape_derivatives(4, i) * node.GetInitialPosition();
    }

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        const Node<3>& node = GetGeometry()[i];
        act.a1 += shape_derivatives(0, i) * node.Coordinates();
        act.a2 += shape_derivatives(1, i) * node.Coordinates();
        act.a11 += shape_derivatives(2, i) * node.Coordinates();
        act.a12 += shape_derivatives(3, i) * node.Coordinates();
        act.a22 += shape_derivatives(4, i) * node.Coordinates();
    }

    Vector3 ref_gab;
    ref_gab(0) = inner_prod(ref.a1, ref.a1);
    ref_gab(1) = inner_prod(ref.a2, ref.a2);
    ref_gab(2) = inner_prod(ref.a1, ref.a2);

    Vector3 act_gab;
    act_gab(0) = inner_prod(act.a1, act.a1);
    act_gab(1) = inner_prod(act.a2, act.a2);
    act_gab(2) = inner_prod(act.a1, act.a2);

    IgaDebug::CheckVector(expected_data, "gab_ref", ref_gab);
    IgaDebug::CheckVector(expected_data, "g1", act.a1);
    IgaDebug::CheckVector(expected_data, "g2", act.a2);

    MathUtils<double>::UnitCrossProduct(ref.a3, ref.a1, ref.a2);
    MathUtils<double>::UnitCrossProduct(act.a3, act.a1, act.a2);

    Vector3 ref_bv;
    ref_bv[0] = inner_prod(ref.a11, ref.a3);
    ref_bv[1] = inner_prod(ref.a12, ref.a3);
    ref_bv[2] = inner_prod(ref.a22, ref.a3);

    Vector3 act_bv;
    act_bv[0] = inner_prod(act.a11, act.a3);
    act_bv[1] = inner_prod(act.a12, act.a3);
    act_bv[2] = inner_prod(act.a22, act.a3);

    IgaDebug::CheckVector(expected_data, "bv_ref", ref_bv);
    IgaDebug::CheckVector(expected_data, "bv", act_bv);

    // transformation matrix

    Vector3 e1 = ref.a1 / MathUtils<double>::Norm(ref.a1);
    Vector3 e2 = ref.a2 - inner_prod(ref.a2, e1) * e1;
    e2 /= MathUtils<double>::Norm(e2);

    Vector3 g_ab_con;
    g_ab_con[0] = ref_gab[1] / (ref_gab[0] * ref_gab[1] - ref_gab[2] * ref_gab[2]);
    g_ab_con[1] = ref_gab[0] / (ref_gab[0] * ref_gab[1] - ref_gab[2] * ref_gab[2]);
    g_ab_con[2] = -ref_gab[2] / (ref_gab[0] * ref_gab[1] - ref_gab[2] * ref_gab[2]);

    Vector3 g_con1 = g_ab_con[0] * ref.a1 + g_ab_con[2] * ref.a2;
    Vector3 g_con2 = g_ab_con[2] * ref.a1 + g_ab_con[1] * ref.a2;

    double eg11 = inner_prod(e1, g_con1);
    double eg12 = inner_prod(e1, g_con2);
    double eg21 = inner_prod(e2, g_con1);
    double eg22 = inner_prod(e2, g_con2);

    BoundedMatrix<double, 3, 3> Tm;
    Tm(0, 0) = eg11 * eg11;
    Tm(0, 1) = eg12 * eg12;
    Tm(0, 2) = 2 * eg11 * eg12;
    Tm(1, 0) = eg21 * eg21;
    Tm(1, 1) = eg22 * eg22;
    Tm(1, 2) = 2 * eg21 * eg22;
    Tm(2, 0) = 2 * eg11 * eg21;
    Tm(2, 1) = 2 * eg12 * eg22;
    Tm(2, 2) = 2 * (eg11 * eg22 + eg12 * eg21);

    IgaDebug::CheckMatrix(expected_data, "Tm", Tm);

    // strain and curvature

    Vector3 E_cu = 0.5 * (act_gab - ref_gab);
    Vector3 K_cu = act_bv - ref_bv;

    Vector3 E_ca = prod(Tm, E_cu);
    Vector3 K_ca = prod(Tm, K_cu);

    IgaDebug::CheckVector(expected_data, "E_cu", E_cu);
    IgaDebug::CheckVector(expected_data, "K_cu", K_cu);
    IgaDebug::CheckVector(expected_data, "E_ca", E_ca);
    IgaDebug::CheckVector(expected_data, "K_ca", K_ca);

    KRATOS_CATCH("")
}

void IgaShell3PElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaShell3PElement\" #" << Id();
}

} // namespace Kratos