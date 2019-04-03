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

    const HyperDualVector<3> tb = b1 / b;
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
    auto const penalty_disp_u = GetValue(PENALTY_DISPLACEMENT_X);
    auto const penalty_disp_v = GetValue(PENALTY_DISPLACEMENT_Y);
    auto const penalty_disp_w = GetValue(PENALTY_DISPLACEMENT_Z);
    auto const penalty_rot    = GetValue(PENALTY_ROTATION);
    auto const penalty_tors   = GetValue(PENALTY_TORSION);

    // const auto dt_1 = ta.dot(a1);
    // const auto dt_2 = b2.dot(ta);
    // const auto dt_3 = b3.dot(ta);

    const double D_11 = Tb.dot(Ta);
    const double D_12 = Tb.dot(A2);
    const double D_13 = Tb.dot(A3);
    const double D_21 = B2.dot(Ta);
    const double D_22 = B2.dot(A2);
    const double D_23 = B2.dot(A3);
    const double D_31 = B3.dot(Ta);
    const double D_32 = B3.dot(A2);
    const double D_33 = B3.dot(A3);

    // std::cout << " D_11 " << D_11 << std::endl;
    // std::cout << " D_12 " << D_12 << std::endl;
    // std::cout << " D_13 " << D_13 << std::endl;


    const auto d_11 = tb.dot(ta);
    const auto d_12 = tb.dot(a2);
    const auto d_13 = tb.dot(a3);
    const auto d_21 = b2.dot(ta);
    const auto d_22 = b2.dot(a2);
    const auto d_23 = b2.dot(a3);
    const auto d_31 = b3.dot(ta);
    const auto d_32 = b3.dot(a2);
    const auto d_33 = b3.dot(a3);

    const auto dP_disp_u = 0.5 * pow(xb_1 - xa_1, 2) * penalty_disp_u;
    const auto dP_disp_v = 0.5 * pow(xb_2 - xa_2, 2) * penalty_disp_v;
    const auto dP_disp_w = 0.5 * pow(xb_3 - xa_3, 2) * penalty_disp_w;

    HyperDualVector<3> e1;
    e1 << D_11 * tb[0] + D_21 * b2[0] + D_31 * b3[0],
          D_11 * tb[1] + D_21 * b2[1] + D_31 * b3[1],
          D_11 * tb[2] + D_21 * b2[2] + D_31 * b3[2];
    HyperDualVector<3> e2;
    e2 << D_12 * tb[0] + D_22 * b2[0] + D_32 * b3[0],
          D_12 * tb[1] + D_22 * b2[1] + D_32 * b3[1],
          D_12 * tb[2] + D_22 * b2[2] + D_32 * b3[2];
    HyperDualVector<3> e3;
    e3 << D_13 * tb[0] + D_23 * b2[0] + D_33 * b3[0],
          D_13 * tb[1] + D_23 * b2[1] + D_33 * b3[1],
          D_13 * tb[2] + D_23 * b2[2] + D_33 * b3[2];

    // const auto alpha_2  = 0.5 * (HyperJet::atan2(d_13, d_11) + HyperJet::atan2(d_31, d_22) // );
    //                             -HyperJet::atan2(D_13, D_11) - HyperJet::atan2(D_31, D_22));   
    // const auto alpha_3  = 0.5 * (HyperJet::atan2(d_12, d_11) + HyperJet::atan2(d_21, d_22)  // );
    //                             -HyperJet::atan2(D_12, D_11) - HyperJet::atan2(D_21, D_22));  
                                
    const auto alpha_112 = HyperJet::atan2(d_11, d_12) - HyperJet::atan2(D_11, D_12); 
    const auto alpha_113 = HyperJet::atan2(d_11, d_13) - HyperJet::atan2(D_11, D_13); 
    const auto alpha_123 = HyperJet::atan2(d_12, d_13) - HyperJet::atan2(D_12, D_13); 
                              
    const auto alpha_212 = HyperJet::atan2(d_21, d_22) - HyperJet::atan2(D_21, D_22); 
    const auto alpha_213 = HyperJet::atan2(d_21, d_23) - HyperJet::atan2(D_21, D_23); 
    const auto alpha_223 = HyperJet::atan2(d_22, d_23) - HyperJet::atan2(D_22, D_23); 
                              
    const auto alpha_312 = HyperJet::atan2(d_31, d_32) - HyperJet::atan2(D_31, D_32); 
    const auto alpha_313 = HyperJet::atan2(d_31, d_33) - HyperJet::atan2(D_31, D_33); 
    const auto alpha_323 = HyperJet::atan2(d_32, d_33) - HyperJet::atan2(D_32, D_33);

    const auto delta_11 = D_11 - d_11;
    const auto delta_12 = D_12 - d_12;
    const auto delta_13 = D_13 - d_13;
    const auto delta_21 = D_21 - d_21;
    const auto delta_22 = D_22 - d_22;
    const auto delta_23 = D_23 - d_23;
    const auto delta_31 = D_31 - d_31;
    const auto delta_32 = D_32 - d_32;
    const auto delta_33 = D_33 - d_33;

    // const auto rotation = 0.5 * (delta_11 * delta_11 + delta_12 * delta_12 + delta_13 * delta_13 + 
    //                              delta_21 * delta_21 + delta_22 * delta_22 + delta_23 * delta_23 + 
    //                              delta_31 * delta_31 + delta_32 * delta_32 + delta_33 * delta_33) * penalty_rot;
    
    const auto d1 = ta - e1;
    const auto d2 = a2.transpose() - e2;
    const auto d3 = a3.transpose() - e3;

    const HyperDualVector<1> rotation = 0.5 * (d1 * d1.transpose() + d2 * d2.transpose() + d3 * d3.transpose()) * penalty_rot;
    // const auto rotation = 0.5 * (alpha_112 * alpha_112 + alpha_113 * alpha_113 + alpha_123 * alpha_123 + 
    //                              alpha_212 * alpha_212 + alpha_213 * alpha_213 + alpha_223 * alpha_223) * penalty_rot;
    // const auto alpha_2  = 0.5 * (HyperJet::atan2(tb.dot(a3), tb.dot(a1)) // );
    //                             -HyperJet::atan2(Tb.dot(A3), Tb.dot(A1)));   
    // const auto alpha_3  = 0.5 * (HyperJet::atan2(tb.dot(a2), tb.dot(a1))  // );
    //                             -HyperJet::atan2(Tb.dot(A2), Tb.dot(A1)));  

    // std::cout << "atan D12 , D11 " << HyperJet::atan2(D_12, D_11) << std::endl;
    // std::cout << "atan D21 , D11 " << HyperJet::atan2(D_21, D_22) << std::endl;




    const auto alpha_12 = (HyperJet::atan2(b2.dot(a3) , b2.dot(a2)) );
                        // -  HyperJet::atan2(B2.dot(A3) , B2.dot(A2)));
    const auto alpha_13 = (HyperJet::atan2(b3.dot(a2) , b3.dot(a3))  );
                        // -  HyperJet::atan2(B3.dot(A2) , B3.dot(A3))); 

    
    
    // // By the Angle
    // const auto dP_alpha_bend_2 = 0.5 * (alpha_2 * alpha_2) * penalty_rot;
    // const auto dP_alpha_bend_3 = 0.5 * (alpha_3 * alpha_3) * penalty_rot;
    // // Potential Torsion
    const auto dP_alpha_tors = 0.5 * 0.5 * (alpha_12 * alpha_12 + alpha_13 * alpha_13) * penalty_tors;
    // const auto dP_alpha_tors = 0.5 * (alpha_12 * alpha_12 ) * penalty_tors;

    const auto f = 0.5 * u.dot(u) * 1e9;

    // MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h() + dP_alpha_bend_2.h() + dP_alpha_bend_3.h() + dP_alpha_tors.h());
    // MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g() + dP_alpha_bend_2.g() + dP_alpha_bend_3.g() + dP_alpha_tors.g());

    
    MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h() + rotation[0].h() );
    MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g() + rotation[0].g() );

    // std::cout << "MapMatrix(rLeftHandSideMatrix) " << std::endl;
    // std::cout <<  MapMatrix(rLeftHandSideMatrix) << std::endl;
    // std::cout << "MapMatrix(rRightHandSideVector) " << std::endl;
    // std::cout <<  MapMatrix(rLeftHandSideMatrix) << std::endl;

    KRATOS_CATCH("")
}

void IgaBeamADWeakCoupling::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaBeamADWeakCoupling\" #" << Id();
}

} // namespace Kratos