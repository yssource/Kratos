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

    Vector3 a1_x_a2;
    MathUtils<double>::CrossProduct(a1_x_a2, ref.a1, ref.a2);
    double dA = MathUtils<double>::Norm(a1_x_a2);

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

    Vector3 dE_cu;
    Vector3 dK_cu;

    vector<double> S_g3dg3(NumberOfDofs());
    vector<double> S_g3dg3lg3_3(NumberOfDofs());
    matrix<double> S_dg_1(3, NumberOfDofs());
    matrix<double> S_dg_2(3, NumberOfDofs());
    matrix<double> S_dg3(3, NumberOfDofs());
    matrix<double> S_dE_ca(3, NumberOfDofs());
    matrix<double> S_dK_ca(3, NumberOfDofs());
    matrix<double> S_dn(3, NumberOfDofs());

    S_g3dg3.clear();
    S_g3dg3lg3_3.clear();
    S_dg_1.clear();
    S_dg_2.clear();
    S_dg3.clear();
    S_dE_ca.clear();
    S_dK_ca.clear();
    S_dn.clear();


    double lg3_3 = dA * dA * dA;
    double lg3_5 = dA * dA * dA * dA * dA;
    double inv_lg3 = 1.0 / dA;
    double inv_lg3_3 = 1.0 / lg3_3;
    double inv_lg3_5 = 1.0 / lg3_5;

    for(std::size_t r = 0; r < NumberOfDofs(); r++) {
        std::size_t shape_index_r = GetShapeIndex(r);
        std::size_t dof_type_index_r = GetDofTypeIndex(r);

        S_dg_1(dof_type_index_r, r) = shape_derivatives(0, shape_index_r);
        S_dg_2(dof_type_index_r, r) = shape_derivatives(1, shape_index_r);

        // strain
        dE_cu[0] = shape_derivatives(0, shape_index_r) * act.a1[dof_type_index_r];
        dE_cu[1] = shape_derivatives(1, shape_index_r) * act.a2[dof_type_index_r];
        dE_cu[2] = 0.5 * (shape_derivatives(0, shape_index_r) * act.a2[dof_type_index_r] + act.a1[dof_type_index_r] * shape_derivatives(1, shape_index_r));

        S_dE_ca(0, r) = Tm(0, 0) * dE_cu[0] + Tm(0, 1) * dE_cu[1] + Tm(0, 2) * dE_cu[2];
        S_dE_ca(1, r) = Tm(1, 0) * dE_cu[0] + Tm(1, 1) * dE_cu[1] + Tm(1, 2) * dE_cu[2];
        S_dE_ca(2, r) = Tm(2, 0) * dE_cu[0] + Tm(2, 1) * dE_cu[1] + Tm(2, 2) * dE_cu[2];

        // curvature
        S_dg3(0, r) = S_dg_1(1, r) * act.a2[2] - S_dg_1(2, r) * act.a2[1] + act.a1(1) * S_dg_2(2, r) - act.a1[2] * S_dg_2(1, r);
        S_dg3(1, r) = S_dg_1(2, r) * act.a2[0] - S_dg_1(0, r) * act.a2[2] + act.a1(2) * S_dg_2(0, r) - act.a1[0] * S_dg_2(2, r);
        S_dg3(2, r) = S_dg_1(0, r) * act.a2[1] - S_dg_1(1, r) * act.a2[0] + act.a1(0) * S_dg_2(1, r) - act.a1[1] * S_dg_2(0, r);

        S_g3dg3[r] = a1_x_a2[0] * S_dg3(0, r)+a1_x_a2[1] * S_dg3(1, r)+a1_x_a2[2] * S_dg3(2, r);
        S_g3dg3lg3_3[r] = S_g3dg3[r] * inv_lg3_3;

        S_dn(0, r) = S_dg3(0, r) * inv_lg3 - a1_x_a2[0] * S_g3dg3lg3_3[r];
        S_dn(1, r) = S_dg3(1, r) * inv_lg3 - a1_x_a2[1] * S_g3dg3lg3_3[r];
        S_dn(2, r) = S_dg3(2, r) * inv_lg3 - a1_x_a2[2] * S_g3dg3lg3_3[r];

        dK_cu[0] = shape_derivatives(2, shape_index_r) * act.a3[dof_type_index_r] + act.a11[0] * S_dn(0, r) + act.a11[1] * S_dn(1, r) + act.a11[2] * S_dn(2, r);
        dK_cu[1] = shape_derivatives(3, shape_index_r) * act.a3[dof_type_index_r] + act.a12[0] * S_dn(0, r) + act.a12[1] * S_dn(1, r) + act.a12[2] * S_dn(2, r);
        dK_cu[2] = shape_derivatives(4, shape_index_r) * act.a3[dof_type_index_r] + act.a22[0] * S_dn(0, r) + act.a22[1] * S_dn(1, r) + act.a22[2] * S_dn(2, r);

        S_dK_ca(0, r) = Tm(0, 0) * dK_cu[0] + Tm(0, 1) * dK_cu[1] + Tm(0, 2) * dK_cu[2];
        S_dK_ca(1, r) = Tm(1, 0) * dK_cu[0] + Tm(1, 1) * dK_cu[1] + Tm(1, 2) * dK_cu[2];
        S_dK_ca(2, r) = Tm(2, 0) * dK_cu[0] + Tm(2, 1) * dK_cu[1] + Tm(2, 2) * dK_cu[2];
    }

    IgaDebug::CheckMatrix(expected_data, "S_dg_1", S_dg_1);
    IgaDebug::CheckMatrix(expected_data, "S_dg_2", S_dg_2);


    matrix<double> S_ddE_ca_1(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddE_ca_2(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddE_ca_3(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddK_ca_1(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddK_ca_2(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddK_ca_3(NumberOfDofs(),NumberOfDofs());
    
    S_ddE_ca_1.clear();
    S_ddE_ca_2.clear();
    S_ddE_ca_3.clear();
    S_ddK_ca_1.clear();
    S_ddK_ca_2.clear();
    S_ddK_ca_3.clear();


    for(std::size_t r = 0; r < NumberOfDofs(); r++) {
        std::size_t shape_index_r = GetShapeIndex(r);
        std::size_t dof_type_index_r = GetDofTypeIndex(r);

        for(std::size_t s = 0; s <= r; s++) {
            std::size_t shape_index_s = GetShapeIndex(s);
            std::size_t dof_type_index_s = GetDofTypeIndex(s);

            // strain
            Vector3 ddE_cu;

            if (dof_type_index_r == dof_type_index_s) {
                ddE_cu[0] = shape_derivatives(0, shape_index_r) * shape_derivatives(0, shape_index_s);
                ddE_cu[1] = shape_derivatives(1, shape_index_r) * shape_derivatives(1, shape_index_s);
                ddE_cu[2] = 0.5 * (shape_derivatives(0, shape_index_r)*shape_derivatives(1, shape_index_s)+shape_derivatives(1, shape_index_r)*shape_derivatives(0, shape_index_s));
            } else {
                ddE_cu.clear();
            }

            S_ddE_ca_1(r, s) = Tm(0, 0) * ddE_cu[0] + Tm(0, 1) * ddE_cu[1] + Tm(0, 2) * ddE_cu[2];
            S_ddE_ca_2(r, s) = Tm(1, 0) * ddE_cu[0] + Tm(1, 1) * ddE_cu[1] + Tm(1, 2) * ddE_cu[2];
            S_ddE_ca_3(r, s) = Tm(2, 0) * ddE_cu[0] + Tm(2, 1) * ddE_cu[1] + Tm(2, 2) * ddE_cu[2];

            // curvature
            Vector3 ddg3;
            ddg3.clear();
            auto dirt = 4-dof_type_index_r-dof_type_index_s;
            auto ddir = dof_type_index_r-dof_type_index_s;
            if(dof_type_index_r + 1 == dof_type_index_s)   ddg3[dirt-1] =  shape_derivatives(0, shape_index_r)*shape_derivatives(1, shape_index_s)-shape_derivatives(0, shape_index_s)*shape_derivatives(1, shape_index_r);
            else if(dof_type_index_r== dof_type_index_s+2) ddg3[dirt-1] =  shape_derivatives(0, shape_index_r)*shape_derivatives(1, shape_index_s)-shape_derivatives(0, shape_index_s)*shape_derivatives(1, shape_index_r);
            else if(dof_type_index_r== dof_type_index_s+1) ddg3[dirt-1] = -shape_derivatives(0, shape_index_r)*shape_derivatives(1, shape_index_s)+shape_derivatives(0, shape_index_s)*shape_derivatives(1, shape_index_r);
            else if(dof_type_index_r + 2 == dof_type_index_s) ddg3[dirt-1] = -shape_derivatives(0, shape_index_r)*shape_derivatives(1, shape_index_s)+shape_derivatives(0, shape_index_s)*shape_derivatives(1, shape_index_r);

            double c = -( ddg3[0]*a1_x_a2[0] + ddg3[1]*a1_x_a2[1] + ddg3[2]*a1_x_a2[2]
                + S_dg3(0,r)*S_dg3(0,s) + S_dg3(1,r)*S_dg3(1,s) + S_dg3(2,r)*S_dg3(2,s)
                )*inv_lg3_3;

            double d = 3.0*S_g3dg3[r]*S_g3dg3[s]*inv_lg3_5;

            Vector3 ddn;
            ddn.clear();
            ddn[0] = ddg3[0] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(0, r) - S_g3dg3lg3_3[r] * S_dg3(0, s) + (c + d) * a1_x_a2[0];
            ddn[1] = ddg3[1] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(1, r) - S_g3dg3lg3_3[r] * S_dg3(1, s) + (c + d) * a1_x_a2[1];
            ddn[2] = ddg3[2] * inv_lg3 - S_g3dg3lg3_3[s] * S_dg3(2, r) - S_g3dg3lg3_3[r] * S_dg3(2, s) + (c + d) * a1_x_a2[2];

            Vector3 ddK_cu;
            ddK_cu.clear();
            ddK_cu[0] = shape_derivatives(2, shape_index_r)*S_dn(dof_type_index_r,s) + shape_derivatives(2, shape_index_s)*S_dn(dof_type_index_s,r)
                        + act.a11[0]*ddn[0] + act.a11[1]*ddn[1] + act.a11[2]*ddn[2];
            ddK_cu[1] = shape_derivatives(3, shape_index_r)*S_dn(dof_type_index_r,s) + shape_derivatives(3, shape_index_s)*S_dn(dof_type_index_s,r)
                        + act.a12[0]*ddn[0] + act.a12[1]*ddn[1] + act.a12[2]*ddn[2];
            ddK_cu[2] = shape_derivatives(4, shape_index_r)*S_dn(dof_type_index_r,s) + shape_derivatives(4, shape_index_s)*S_dn(dof_type_index_s,r)
                        + act.a22[0]*ddn[0] + act.a22[1]*ddn[1] + act.a22[2]*ddn[2];

            S_ddK_ca_1(r,s) = Tm(0,0)*ddK_cu[0] + Tm(0,1)*ddK_cu[1] + Tm(0,2)*ddK_cu[2];
            S_ddK_ca_2(r,s) = Tm(1,0)*ddK_cu[0] + Tm(1,1)*ddK_cu[1] + Tm(1,2)*ddK_cu[2];
            S_ddK_ca_3(r,s) = Tm(2,0)*ddK_cu[0] + Tm(2,1)*ddK_cu[1] + Tm(2,2)*ddK_cu[2];
        }
    }

    IgaDebug::CheckLowerMatrix(expected_data, "S_ddE_ca_1", S_ddE_ca_1);
    IgaDebug::CheckLowerMatrix(expected_data, "S_ddE_ca_2", S_ddE_ca_2);
    IgaDebug::CheckLowerMatrix(expected_data, "S_ddE_ca_3", S_ddE_ca_3);

    IgaDebug::CheckLowerMatrix(expected_data, "S_ddK_ca_1", S_ddK_ca_1);
    IgaDebug::CheckLowerMatrix(expected_data, "S_ddK_ca_2", S_ddK_ca_2);
    IgaDebug::CheckLowerMatrix(expected_data, "S_ddK_ca_3", S_ddK_ca_3);


    matrix<double> S_kem(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_keb(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_dK(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddK(NumberOfDofs(),NumberOfDofs());
    matrix<double> S_ddK_1(NumberOfDofs(),NumberOfDofs());
    vector<double> S_fiem(NumberOfDofs());
    vector<double> S_fieb(NumberOfDofs());

    Vector3 N_ca;
    Vector3 M_ca;

    axpy_prod(dm, E_ca, N_ca, true);
    axpy_prod(db, K_ca, M_ca, true);

    for(std::size_t r = 0; r < NumberOfDofs(); r++) {
        Vector3 dN_ca;
        dN_ca[0] = dm(0,0)*S_dE_ca(0,r) + dm(0,1)*S_dE_ca(1,r) + dm(0,2)*S_dE_ca(2,r);
        dN_ca[1] = dm(1,0)*S_dE_ca(0,r) + dm(1,1)*S_dE_ca(1,r) + dm(1,2)*S_dE_ca(2,r);
        dN_ca[2] = dm(2,0)*S_dE_ca(0,r) + dm(2,1)*S_dE_ca(1,r) + dm(2,2)*S_dE_ca(2,r);

        Vector3 dM_ca;
        dM_ca[0] = db(0,0)*S_dK_ca(0,r) + db(0,1)*S_dK_ca(1,r) + db(0,2)*S_dK_ca(2,r);
        dM_ca[1] = db(1,0)*S_dK_ca(0,r) + db(1,1)*S_dK_ca(1,r) + db(1,2)*S_dK_ca(2,r);
        dM_ca[2] = db(2,0)*S_dK_ca(0,r) + db(2,1)*S_dK_ca(1,r) + db(2,2)*S_dK_ca(2,r);

        for(std::size_t s = 0; s < NumberOfDofs(); s++) {
            // membrane stiffness
            S_kem(r,s) = dN_ca[0]*S_dE_ca(0,s) + dN_ca[1]*S_dE_ca(1,s) + dN_ca[2]*S_dE_ca(2,s)
                    + N_ca[0]*S_ddE_ca_1(r,s) + N_ca[1]*S_ddE_ca_2(r,s) + N_ca[2]*S_ddE_ca_3(r,s);
            S_kem(s,r) = S_kem(r,s);
            // bending stiffness
            S_keb(r,s) = dM_ca[0]*S_dK_ca(0,s) + dM_ca[1]*S_dK_ca(1,s) + dM_ca[2]*S_dK_ca(2,s)
                    + M_ca[0]*S_ddK_ca_1(r,s) + M_ca[1]*S_ddK_ca_2(r,s) + M_ca[2]*S_ddK_ca_3(r,s);
            S_keb(s,r) = S_keb(r,s);

            S_dK(r,s)=dM_ca[0]*S_dK_ca(0,s) + dM_ca[1]*S_dK_ca(1,s) + dM_ca[2]*S_dK_ca(2,s);
            S_dK(s,r) = S_dK(r,s);

            S_ddK(r,s)=M_ca[0]*S_ddK_ca_1(r,s) + M_ca[1]*S_ddK_ca_2(r,s) + M_ca[2]*S_ddK_ca_3(r,s);
            S_ddK(s,r) = S_ddK(r,s);

            S_ddK_1(r,s)=M_ca[0]*S_ddK_ca_1(r,s);
            S_ddK_1(s,r) = S_ddK_1(r,s);
            S_ddK_ca_1(s,r) = S_ddK_ca_1(r,s);
            S_ddK_ca_2(s,r) = S_ddK_ca_2(r,s);
            S_ddK_ca_3(s,r) = S_ddK_ca_3(r,s);
        }

        S_fiem(r) = -(N_ca[0]*S_dE_ca(0,r) + N_ca[1]*S_dE_ca(1,r) + N_ca[2]*S_dE_ca(2,r));
        S_fieb(r) = -(M_ca[0]*S_dK_ca(0,r) + M_ca[1]*S_dK_ca(1,r) + M_ca[2]*S_dK_ca(2,r));
   }

    matrix<double> gke = (S_kem+S_keb)*dA;
    vector<double> fie = (S_fiem+S_fieb)*dA;

    IgaDebug::CheckLowerMatrix(expected_data, "_gke", gke);
    IgaDebug::CheckVector(expected_data, "_fie", fie);

    KRATOS_CATCH("")
}

void IgaShell3PElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaShell3PElement\" #" << Id();
}

} // namespace Kratos