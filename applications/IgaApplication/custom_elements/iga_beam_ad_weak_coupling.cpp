/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
//                  Lukas Rauch
*/

// System includes
#include "includes/define.h"
#include "includes/variables.h"
// External includes

// Project includes
#include "iga_beam_ad_weak_coupling.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaBeamADWeakCoupling::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaBeamADWeakCoupling>(NewId, geometry,
        pProperties);
}

void IgaBeamADWeakCoupling::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
        SetElementDof(rElementalDofList, i, 3, DISPLACEMENT_ROTATION);
    }

    KRATOS_CATCH("")
}

void IgaBeamADWeakCoupling::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetElementEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetElementEquationId(rResult, i, 2, DISPLACEMENT_Z);
        SetElementEquationId(rResult, i, 3, DISPLACEMENT_ROTATION);
    }

    KRATOS_CATCH("")
}

void IgaBeamADWeakCoupling::Initialize()
{

}

void IgaBeamADWeakCoupling::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    using namespace BeamUtilities;
    using std::sqrt;
    using std::pow;

    using Vector3d = BeamUtilities::Vector<3>;
    using Matrix3d = BeamUtilities::Matrix<3, 3>;

    KRATOS_TRY;

    // get integration weight

    // const double& integration_weight = GetValue(INTEGRATION_WEIGHT);

    // get shape functions

    BeamUtilities::Matrix<3, Dynamic> shape_functions_a(3, GetValue(SHAPE_FUNCTION_VALUES).size());
    shape_functions_a.row(0) = MapVector(GetValue(SHAPE_FUNCTION_VALUES));
    shape_functions_a.row(1) = MapVector(GetValue(SHAPE_FUNCTION_LOCAL_DER_1));
    shape_functions_a.row(2) = MapVector(GetValue(SHAPE_FUNCTION_LOCAL_DER_2));

    const Vector3d A1   = MapVector(GetValue(BASE_A1));
    const Vector3d A1_1 = MapVector(GetValue(BASE_A1_1));

    const double A11 = A1.dot(A1);
    const double A = sqrt(A11);

    const Vector3d Ta = A1 / A;
    const Vector3d Ta_1 = A1_1 / A - A1.dot(A1_1) * A1 / pow(A, 3);

    const Vector3d A2   = MapVector(GetValue(BASE_A2));
    const Vector3d A3   = MapVector(GetValue(BASE_A3));

    const auto phia = ComputeActValue(DISPLACEMENT_ROTATION, 0, shape_functions_a, GetGeometry());
    const auto phia_1 = ComputeActValue(DISPLACEMENT_ROTATION, 1, shape_functions_a, GetGeometry());    

    const auto xa = ComputeActBaseVector(0, shape_functions_a, GetGeometry());
    const auto a1 = ComputeActBaseVector(1, shape_functions_a, GetGeometry());
    const auto a1_1 = ComputeActBaseVector(2, shape_functions_a, GetGeometry());

    const auto a11 = a1.dot(a1);
    const auto a = sqrt(a11);

    const auto ta = a1 / a;
    const auto ta_1 = a1_1 / a - a1.dot(a1_1) * a1 / pow(a, 3);

    const auto roda = ComputeRod<HyperDual>(ta, phia);
    const auto roda_1 = ComputeRod_1<HyperDual>(ta, ta_1, phia, phia_1);

    const auto lama = ComputeLam<HyperDual>(Ta, ta);
    const auto lama_1 = ComputeLam_1<HyperDual>(Ta, Ta_1, ta, ta_1);

    const auto rod_lama = roda * lama;

    const auto a2 = rod_lama * A2.transpose();
    const auto a3 = rod_lama * A3.transpose();


    BeamUtilities::Matrix<3, Dynamic> shape_functions_b(3, GetValue(SHAPE_FUNCTION_VALUES_B).size());
    shape_functions_b.row(0) = MapVector(GetValue(SHAPE_FUNCTION_VALUES_B));
    shape_functions_b.row(1) = MapVector(GetValue(SHAPE_FUNCTION_LOCAL_DER_1_B));
    shape_functions_b.row(2) = MapVector(GetValue(SHAPE_FUNCTION_LOCAL_DER_2_B));

    const Vector3d B1   = MapVector(GetValue(BASE_B1));
    const Vector3d B1_1 = MapVector(GetValue(BASE_B1_1));

    const double B11 = B1.dot(B1);
    const double B = sqrt(B11);

    const Vector3d Tb = B1 / B;
    const Vector3d Tb_1 = B1_1 / B - B1.dot(B1_1) * B1 / pow(B, 3);

    const Vector3d B2   = MapVector(GetValue(BASE_B2));
    const Vector3d B3   = MapVector(GetValue(BASE_B3));

    const auto phib   = ComputeActValue(DISPLACEMENT_ROTATION, 0, shape_functions_b, GetGeometry(), true);
    const auto phib_1 = ComputeActValue(DISPLACEMENT_ROTATION, 1, shape_functions_b, GetGeometry(), true);    

    const auto xb = ComputeActBaseVector(0, shape_functions_b, GetGeometry(), true);
    const auto b1 = ComputeActBaseVector(1, shape_functions_b, GetGeometry(), true);
    const auto b1_1 = ComputeActBaseVector(2, shape_functions_b, GetGeometry(), true);

    const auto b11 = b1.dot(b1);
    const auto b = sqrt(b11);

    const auto tb = b1 / b;
    const auto tb_1 = b1_1 / b - b1.dot(b1_1) * b1 / pow(b, 3);

    const auto rodb = ComputeRod<HyperDual>(tb, phib);
    const auto rodb_1 = ComputeRod_1<HyperDual>(tb, tb_1, phib, phib_1);

    const auto lamb = ComputeLam<HyperDual>(Tb, tb);
    const auto lamb_1 = ComputeLam_1<HyperDual>(Tb, Tb_1, tb, tb_1);

    const auto rod_lamb = rodb * lamb;

    const auto b2 = rod_lamb * B2.transpose();
    const auto b3 = rod_lamb * B3.transpose();


    const auto xa_1 = xa.dot(ta);
    const auto xa_2 = xa.dot(a2);
    const auto xa_3 = xa.dot(a3);

    const auto xb_1 = xb.dot(ta);
    const auto xb_2 = xb.dot(a2);
    const auto xb_3 = xb.dot(a3);


    const auto u = xb - xa;

    
    // // penalties
    auto const penalty_disp_u = 1e9; //GetValue(PENALTY_DISPLACEMENT_X);
    auto const penalty_disp_v = 1e9; //GetValue(PENALTY_DISPLACEMENT_Y);
    auto const penalty_disp_w = 1e9; //GetValue(PENALTY_DISPLACEMENT_Z);
    auto const penalty_rot    = 1e9; //GetValue(PENALTY_ROTATION);
    auto const penalty_tors   = 1e9; //GetValue(PENALTY_TORSION);

    const auto dt_1 = ta.dot(a1);
    const auto dt_2 = b2.dot(ta);
    const auto dt_3 = b3.dot(ta);


    const auto d_11 = tb.dot(ta);
    const auto d_12 = tb.dot(a2);
    const auto d_13 = tb.dot(a3);
    const auto d_21 = b2.dot(a1);
    const auto d_22 = b2.dot(a2);
    const auto d_31 = b3.dot(a1);
    const auto d_33 = b3.dot(a3);

    const auto d_t = tb.dot(ta);     // Anteil T in a1 Richtung
    const auto d_n = tb.dot(a2);     // Anteil N in a1 Richtung
    const auto d_v = tb.dot(a3);     // Anteil V in a1 Richtung

    // const auto alpha_12 = HyperJet::atan2(a2.dot(A3) , a2.dot(A2));    // Winkel zwischen a2 und A
    // const auto alpha_13 = HyperJet::atan2(a3.dot(A2) , a3.dot(A3));    // Winkel zwischen a2 und A
    // const auto alpha_2  = HyperJet::atan2(d_n , d_t);                  // Winkel zwischen a1 und A1 um n
    // const auto alpha_3  = HyperJet::atan2(d_v , d_t);                  // Winkel zwischen a1 und A1 um v


    // const auto alpha_12 = HyperJet::atan2(b2.dot(a3) , b2.dot(a2));    
    // const auto alpha_13 = HyperJet::atan2(b3.dot(a2) , b3.dot(a3));    

    const auto alpha_2  = 0.5 * (HyperJet::atan2(d_13, d_11) + HyperJet::atan2(d_31, d_22));   
    const auto alpha_3  = 0.5 * (HyperJet::atan2(d_12, d_11) + HyperJet::atan2(d_21, d_22));   

    const auto dP_disp_u = 0.5 * pow(xb_1 - xa_1, 2) * penalty_disp_u;
    const auto dP_disp_v = 0.5 * pow(xb_2 - xa_2, 2) * penalty_disp_v;
    const auto dP_disp_w = 0.5 * pow(xb_3 - xa_3, 2) * penalty_disp_w;
  
    // const auto dP_disp_u = 0.5 * pow(xb(0) - xa(0), 2) * penalty_disp_u;
    // const auto dP_disp_v = 0.5 * pow(xb(1) - xa(1), 2) * penalty_disp_v;
    // const auto dP_disp_w = 0.5 * pow(xb(2) - xa(2), 2) * penalty_disp_w;

    // // By the Angle
    const auto dP_alpha_bend_2 = 0.5 * (alpha_2 * alpha_2) * penalty_rot;
    const auto dP_alpha_bend_3 = 0.5 * (alpha_3 * alpha_3) * penalty_rot;
    // // Potential Torsion
    // const auto dP_alpha_tors = 0.5 * (alpha_12 * alpha_12 + alpha_13 * alpha_13) * penalty_tors;

    // By the Angle
    // const auto dP_alpha_bend = 0.5 * (pow(HyperJet::acos(tb.dot(ta)), 2)) * penalty_rot;
    // Potential Torsion
    // const auto dP_alpha_tors = 0.5 * (alpha_12 * alpha_12 + alpha_13 * alpha_13) * penalty_tors;

    const auto f = 0.5 * u.dot(u) * 1e9;

    // std::cout << "displac u: " <<  dP_disp_u.g() << std::endl; 
    // std::cout << "grad bend: " <<  dP_alpha_bend.g() << std::endl; 
    // std::cout << "hess bend: " <<  dP_alpha_bend.h() << std::endl; 

    // MapMatrix(rLeftHandSideMatrix) = f.h(); //( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h()+ dP_alpha_tors.h() + dP_alpha_bend.h()) ;
    // MapVector(rRightHandSideVector) = -f.g();//-( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g()+ dP_alpha_tors.g() + dP_alpha_bend.g()) ;
    MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h() + dP_alpha_bend_2.h() + dP_alpha_bend_3.h());
    MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g() + dP_alpha_bend_2.g() + dP_alpha_bend_3.g());
    // MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h()+ dP_alpha_tors.h() + dP_alpha_bend.h()) ;
    // MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g()+ dP_alpha_tors.g() + dP_alpha_bend.g()) ;

    KRATOS_CATCH("")
}

void IgaBeamADWeakCoupling::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaBeamADWeakCoupling\" #" << Id();
}

} // namespace Kratos