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
#include "custom_utilities/grid.h"
#include "custom_utilities/iga_debug.h"

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

    m_reference_configuration = ComputeReferenceConfiguration();

    KRATOS_CATCH("")
}

IgaShell3PElement::Configuration
IgaShell3PElement::ComputeReferenceConfiguration()
{
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    Configuration ref;

    ref.a1.clear();
    ref.a2.clear();
    ref.a3.clear();
    ref.a11.clear();
    ref.a12.clear();
    ref.a22.clear();

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        const Node<3>& node = GetGeometry()[i];

        ref.a1 += shape_derivatives(0, i) * node.GetInitialPosition();
        ref.a2 += shape_derivatives(1, i) * node.GetInitialPosition();
        ref.a11 += shape_derivatives(2, i) * node.GetInitialPosition();
        ref.a12 += shape_derivatives(3, i) * node.GetInitialPosition();
        ref.a22 += shape_derivatives(4, i) * node.GetInitialPosition();
    }

    MathUtils<double>::UnitCrossProduct(ref.a3, ref.a1, ref.a2);

    ref.gab[0] = inner_prod(ref.a1, ref.a1);
    ref.gab[1] = inner_prod(ref.a2, ref.a2);
    ref.gab[2] = inner_prod(ref.a1, ref.a2);

    ref.bv[0] = inner_prod(ref.a11, ref.a3);
    ref.bv[1] = inner_prod(ref.a12, ref.a3);
    ref.bv[2] = inner_prod(ref.a22, ref.a3);

    return ref;
}

IgaShell3PElement::Configuration
IgaShell3PElement::ComputeActualConfiguration()
{
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    Configuration act;

    act.a1.clear();
    act.a2.clear();
    act.a3.clear();
    act.a11.clear();
    act.a12.clear();
    act.a22.clear();

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        const Node<3>& node = GetGeometry()[i];

        act.a1 += shape_derivatives(0, i) * node.Coordinates();
        act.a2 += shape_derivatives(1, i) * node.Coordinates();
        act.a11 += shape_derivatives(2, i) * node.Coordinates();
        act.a12 += shape_derivatives(3, i) * node.Coordinates();
        act.a22 += shape_derivatives(4, i) * node.Coordinates();
    }

    MathUtils<double>::UnitCrossProduct(act.a3, act.a1, act.a2);

    act.gab[0] = inner_prod(act.a1, act.a1);
    act.gab[1] = inner_prod(act.a2, act.a2);
    act.gab[2] = inner_prod(act.a1, act.a2);

    act.bv[0] = inner_prod(act.a11, act.a3);
    act.bv[1] = inner_prod(act.a12, act.a3);
    act.bv[2] = inner_prod(act.a22, act.a3);

    return act;
}

void IgaShell3PElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

    #ifdef KRATOS_DEBUG
    // temporary debug data
    auto expected_data = Parameters(GetValue(DEBUG_EXPECTED_DATA));
    #endif

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
    {
        db.clear();
        db(0, 0) = 1.0;
        db(0, 1) = poisson_ratio;
        db(1, 0) = poisson_ratio;
        db(1, 1) = 1.0;
        db(2, 2) = (1.0 - poisson_ratio) / 2.0;
        db *= E * pow(thickness, 3) /
            (12.0 * (1.0 - poisson_ratio * poisson_ratio));
    }

    BoundedMatrix<double, 3, 3> dm;
    {
        dm.clear();
        dm(0, 0) = 1.0;
        dm(0, 1) = poisson_ratio;
        dm(1, 0) = poisson_ratio;
        dm(1, 1) = 1.0;
        dm(2, 2) = (1.0 - poisson_ratio) / 2.0;
        dm *= E * thickness / (1.0 - poisson_ratio * poisson_ratio);
    }

    #ifdef KRATOS_DEBUG
    IgaDebug::CheckMatrix(expected_data, "_Db", db);
    IgaDebug::CheckMatrix(expected_data, "_Dm", dm);
    #endif

    // reference and actual configuration

    Configuration& ref = m_reference_configuration;
    Configuration act = ComputeActualConfiguration();

    #ifdef KRATOS_DEBUG
    IgaDebug::CheckVector(expected_data, "gab_ref", ref.gab);
    IgaDebug::CheckVector(expected_data, "g1", act.a1);
    IgaDebug::CheckVector(expected_data, "g2", act.a2);

    IgaDebug::CheckVector(expected_data, "bv_ref", ref.bv);
    IgaDebug::CheckVector(expected_data, "bv", act.bv);
    #endif

    Vector3 a1_x_a2 = Cross(ref.a1, ref.a2);
    double dA = MathUtils<double>::Norm(a1_x_a2);

    // transformation matrix

    BoundedMatrix<double, 3, 3> Tm;
    {
        Vector3 e1 = ref.a1 / MathUtils<double>::Norm(ref.a1);
        Vector3 e2 = ref.a2 - inner_prod(ref.a2, e1) * e1;
        e2 /= MathUtils<double>::Norm(e2);

        const double det = ref.gab[0] * ref.gab[1] - ref.gab[2] * ref.gab[2];

        Vector3 g_ab_con;
        g_ab_con[0] = ref.gab[1] / det;
        g_ab_con[1] = ref.gab[0] / det;
        g_ab_con[2] = -ref.gab[2] / det;

        const Vector3 g_con1 = g_ab_con[0] * ref.a1 + g_ab_con[2] * ref.a2;
        const Vector3 g_con2 = g_ab_con[2] * ref.a1 + g_ab_con[1] * ref.a2;

        const double eg11 = inner_prod(e1, g_con1);
        const double eg12 = inner_prod(e1, g_con2);
        const double eg21 = inner_prod(e2, g_con1);
        const double eg22 = inner_prod(e2, g_con2);

        Tm(0, 0) = eg11 * eg11;
        Tm(0, 1) = eg12 * eg12;
        Tm(0, 2) = 2 * eg11 * eg12;
        Tm(1, 0) = eg21 * eg21;
        Tm(1, 1) = eg22 * eg22;
        Tm(1, 2) = 2 * eg21 * eg22;
        Tm(2, 0) = 2 * eg11 * eg21;
        Tm(2, 1) = 2 * eg12 * eg22;
        Tm(2, 2) = 2 * (eg11 * eg22 + eg12 * eg21);
    }

    #ifdef KRATOS_DEBUG
    IgaDebug::CheckMatrix(expected_data, "Tm", Tm);
    #endif

    // strain and curvature

    const Vector3 E_cu = 0.5 * (act.gab - ref.gab);
    const Vector3 K_cu = act.bv - ref.bv;

    const Vector3 E_ca = prod(Tm, E_cu);
    const Vector3 K_ca = prod(Tm, K_cu);

    #ifdef KRATOS_DEBUG
    IgaDebug::CheckVector(expected_data, "E_cu", E_cu);
    IgaDebug::CheckVector(expected_data, "K_cu", K_cu);
    IgaDebug::CheckVector(expected_data, "E_ca", E_ca);
    IgaDebug::CheckVector(expected_data, "K_ca", K_ca);
    #endif

    vector<double> S_g3dg3(NumberOfDofs());
    vector<double> S_g3dg3lg3_3(NumberOfDofs());
    std::vector<Vector3> S_dg3(NumberOfDofs());
    std::vector<Vector3> S_dE_ca(NumberOfDofs());
    std::vector<Vector3> S_dK_ca(NumberOfDofs());
    std::vector<Vector3> S_dn(NumberOfDofs());

    Grid<Vector3> S_ddE_ca(NumberOfDofs(), NumberOfDofs());
    Grid<Vector3> S_ddK_ca(NumberOfDofs(), NumberOfDofs());

    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        const std::size_t shape_index_r = GetShapeIndex(r);
        const std::size_t dof_type_index_r = GetDofTypeIndex(r);

        Vector3 S_dg_1;
        Vector3 S_dg_2;

        S_dg_1.clear();
        S_dg_2.clear();

        S_dg_1[dof_type_index_r] = shape_derivatives(0, shape_index_r);
        S_dg_2[dof_type_index_r] = shape_derivatives(1, shape_index_r);

        // strain
        Vector3 dE_cu;
        dE_cu[0] = shape_derivatives(0, shape_index_r) * act.a1[dof_type_index_r];
        dE_cu[1] = shape_derivatives(1, shape_index_r) * act.a2[dof_type_index_r];
        dE_cu[2] = 0.5 * (
                   shape_derivatives(0, shape_index_r) * act.a2[dof_type_index_r] +
                   shape_derivatives(1, shape_index_r) * act.a1[dof_type_index_r]);

        S_dE_ca[r] = prod(Tm, dE_cu);

        // curvature
        S_dg3[r] = Cross(S_dg_1, act.a2) + Cross(act.a1, S_dg_2);

        S_g3dg3[r] = inner_prod(a1_x_a2, S_dg3[r]);
        S_g3dg3lg3_3[r] = S_g3dg3[r] / std::pow(dA, 3);

        S_dn[r] = S_dg3[r] / dA - a1_x_a2 * S_g3dg3lg3_3[r];

        Vector3 dK_cu;
        dK_cu[0] = shape_derivatives(2, shape_index_r) * act.a3[dof_type_index_r] +
                   inner_prod(act.a11, S_dn[r]);
        dK_cu[1] = shape_derivatives(3, shape_index_r) * act.a3[dof_type_index_r] +
                   inner_prod(act.a12, S_dn[r]);
        dK_cu[2] = shape_derivatives(4, shape_index_r) * act.a3[dof_type_index_r] +
                   inner_prod(act.a22, S_dn[r]);

        S_dK_ca[r] = prod(Tm, dK_cu);
    }

    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        const std::size_t shape_index_r = GetShapeIndex(r);
        const std::size_t dof_type_index_r = GetDofTypeIndex(r);

        for (std::size_t s = 0; s <= r; s++) {
            const std::size_t shape_index_s = GetShapeIndex(s);
            const std::size_t dof_type_index_s = GetDofTypeIndex(s);

            // strain

            if (dof_type_index_r == dof_type_index_s) {
                Vector3 ddE_cu;

                ddE_cu[0] = shape_derivatives(0, shape_index_r) *
                            shape_derivatives(0, shape_index_s);
                ddE_cu[1] = shape_derivatives(1, shape_index_r) *
                            shape_derivatives(1, shape_index_s);
                ddE_cu[2] = 0.5 * (
                            shape_derivatives(0, shape_index_r) *
                            shape_derivatives(1, shape_index_s) +
                            shape_derivatives(1, shape_index_r) *
                            shape_derivatives(0, shape_index_s));

                S_ddE_ca(r, s) = prod(Tm, ddE_cu);
            } else {
                S_ddE_ca(r, s).clear();
            }

            // curvature

            Vector3 ddg3;
            ddg3.clear();

            if (dof_type_index_r != dof_type_index_s) {
                const std::size_t ddg3_i =
                    3 - dof_type_index_r - dof_type_index_s;

                const double ddg3_value =
                    shape_derivatives(0, shape_index_r) *
                    shape_derivatives(1, shape_index_s) -
                    shape_derivatives(0, shape_index_s) *
                    shape_derivatives(1, shape_index_r);

                if ((dof_type_index_s == dof_type_index_r + 1) ||
                    (dof_type_index_r == dof_type_index_s + 2)) {
                    ddg3[ddg3_i] = ddg3_value;
                } else {
                    ddg3[ddg3_i] = -ddg3_value;
                }
            }

            const double c = -(inner_prod(ddg3, a1_x_a2) +
                               inner_prod(S_dg3[r], S_dg3[s])) / std::pow(dA, 3);

            const double d = 3.0 * S_g3dg3[r] * S_g3dg3[s] / std::pow(dA, 5);

            const Vector3 ddn = ddg3 / dA - S_g3dg3lg3_3[s] * S_dg3[r] -
                S_g3dg3lg3_3[r] * S_dg3[s] + (c + d) * a1_x_a2;

            Vector3 ddK_cu;
            ddK_cu[0] = shape_derivatives(2, shape_index_r) *
                S_dn[s][dof_type_index_r] +
                shape_derivatives(2, shape_index_s) *
                S_dn[r][dof_type_index_s] +
                act.a11[0] * ddn[0] +
                act.a11[1] * ddn[1] +
                act.a11[2] * ddn[2];
            ddK_cu[1] = shape_derivatives(3, shape_index_r) *
                S_dn[s][dof_type_index_r] +
                shape_derivatives(3, shape_index_s) *
                S_dn[r][dof_type_index_s] +
                act.a12[0] * ddn[0] +
                act.a12[1] * ddn[1] +
                act.a12[2] * ddn[2];
            ddK_cu[2] = shape_derivatives(4, shape_index_r) *
                S_dn[s][dof_type_index_r] +
                shape_derivatives(4, shape_index_s) *
                S_dn[r][dof_type_index_s] +
                act.a22[0] * ddn[0] +
                act.a22[1] * ddn[1] +
                act.a22[2] * ddn[2];

            S_ddK_ca(r, s) = prod(Tm, ddK_cu);
        }
    }

    #ifdef KRATOS_DEBUG
    IgaDebug::CheckLowerGridComponent(expected_data, "S_ddE_ca_1", S_ddE_ca, 0);
    IgaDebug::CheckLowerGridComponent(expected_data, "S_ddE_ca_2", S_ddE_ca, 1);
    IgaDebug::CheckLowerGridComponent(expected_data, "S_ddE_ca_3", S_ddE_ca, 2);

    IgaDebug::CheckLowerGridComponent(expected_data, "S_ddK_ca_1", S_ddK_ca, 0);
    IgaDebug::CheckLowerGridComponent(expected_data, "S_ddK_ca_2", S_ddK_ca, 1);
    IgaDebug::CheckLowerGridComponent(expected_data, "S_ddK_ca_3", S_ddK_ca, 2);
    #endif

    const Vector3 N_ca = prod(dm, E_ca);
    const Vector3 M_ca = prod(db, K_ca);

    for(std::size_t r = 0; r < NumberOfDofs(); r++) {
        const Vector3 dN_ca = prod(dm, S_dE_ca[r]);
        const Vector3 dM_ca = prod(db, S_dK_ca[r]);

        if (ComputeLeftHandSide) {
            for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                // membrane stiffness
                const double S_kem = inner_prod(dN_ca, S_dE_ca[s]) +
                                     inner_prod(N_ca, S_ddE_ca(r, s));

                // bending stiffness
                const double S_keb = inner_prod(dM_ca, S_dK_ca[s]) +
                                     inner_prod(M_ca, S_ddK_ca(r, s));

                rLeftHandSideMatrix(r, s) = (S_kem + S_keb) * dA;

                // symmetry
                rLeftHandSideMatrix(s, r) = rLeftHandSideMatrix(r, s);
            }
        }

        if (ComputeRightHandSide) {
            rRightHandSideVector[r] = -dA * (inner_prod(N_ca, S_dE_ca[r]) +
                inner_prod(M_ca, S_dK_ca[r]));
        }
    }

    #ifdef KRATOS_DEBUG
    if (ComputeLeftHandSide) {
        IgaDebug::CheckLowerMatrix(expected_data, "_gke", rLeftHandSideMatrix);
    }

    if (ComputeRightHandSide) {
        IgaDebug::CheckVector(expected_data, "_fie", rRightHandSideVector);
    }
    #endif

    KRATOS_CATCH("")
}

void IgaShell3PElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaShell3PElement\" #" << Id();
}

} // namespace Kratos