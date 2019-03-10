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
#include "iga_beam_weak_dirichlet_condition.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaBeamWeakDirichletCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaBeamWeakDirichletCondition>(NewId, geometry,
        pProperties);
}

void IgaBeamWeakDirichletCondition::GetDofList(
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

void IgaBeamWeakDirichletCondition::EquationIdVector(
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

void IgaBeamWeakDirichletCondition::Initialize()
{

}

void IgaBeamWeakDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    using namespace BeamUtilities;
    using std::sqrt;
    using std::pow;
    using std::atan;

    using Vector3d = BeamUtilities::Vector<3>;
    using Matrix3d = BeamUtilities::Matrix<3, 3>;

    KRATOS_TRY;

    // get integration weight

    const double& integration_weight = GetValue(INTEGRATION_WEIGHT);

    // get shape functions

    BeamUtilities::Matrix<3, Dynamic> shape_functions(3, NumberOfNodes());
    shape_functions.row(0) = MapVector(GetValue(SHAPE_FUNCTION_VALUES));
    shape_functions.row(1) = MapVector(GetValue(SHAPE_FUNCTION_LOCAL_DER_1));
    shape_functions.row(2) = MapVector(GetValue(SHAPE_FUNCTION_LOCAL_DER_2));
    // shape_functions.bottomRows<2>() = MapMatrix(GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));

    // get properties

    const auto& properties = GetProperties();

    const double young_modulus = properties[YOUNG_MODULUS];
    const double shear_modulus = properties[SHEAR_MODULUS];
    const double area = properties[CROSS_AREA];
    const double moment_of_inertia_x = properties[MOMENT_OF_INERTIA_T];
    const double moment_of_inertia_y = properties[MOMENT_OF_INERTIA_Y];
    const double moment_of_inertia_z = properties[MOMENT_OF_INERTIA_Z];
    // const double prestress = properties[PRESTRESS_CAUCHY];
    // const double Phi = GetValue(PHI);
    // const double Phi_1 = GetValue(PHI_DER_1);

    // const Vector3d A01 = MapVector(GetValue(T0));
    // const Vector3d A01_1 = MapVector(GetValue(T0_DER));
    // const Vector3d A02 = MapVector(GetValue(N0));
    // const Vector3d A03 = A01.cross(A02);

    // material

    const double ea  = young_modulus * area;
    const double gi1 = shear_modulus * moment_of_inertia_x;
    const double ei2 = young_modulus * moment_of_inertia_y;
    const double ei3 = young_modulus * moment_of_inertia_z;

    // reference configuration FIXME: move this section to Initialize()
    const auto X = ComputeRefBaseVector(0, shape_functions, GetGeometry());
    // const Vector3d A1 = ComputeRefBaseVector(1, shape_functions, GetGeometry());
    // const Vector3d A1_1 = ComputeRefBaseVector(2, shape_functions, GetGeometry());
    const Vector3d A1   = MapVector(GetValue(BASE_A1));
    // const Vector3d A1_1 = MapVector(GetValue(BASE_A1_1));

    const double A11 = A1.dot(A1);
    const double A = sqrt(A11);

    const Vector3d T = A1 / A;
    // const Vector3d T_1 = A1_1 / A - A1.dot(A1_1) * A1 / pow(A, 3);

    // const Matrix3d Rod = ComputeRod<double>(T, Phi);
    // const Matrix3d Rod_1 = ComputeRod_1<double>(T, T_1, Phi, Phi_1);

    // const Matrix3d Lam = ComputeLam<double>(A01, T);
    // const Matrix3d Lam_1 = ComputeLam_1<double>(A01, A01_1, T, T_1);

    // const Matrix3d Rod_Lam = Rod * Lam;
    // const Matrix3d Rod_1_Lam = Rod_1 * Lam;
    // const Matrix3d Rod_Lam_1 = Rod * Lam_1;

    // const Vector3d A2 = Rod_Lam * A02.transpose();
    // const Vector3d A2_1 = Rod_1_Lam * A02.transpose() + Rod_Lam_1 * A02.transpose();
    const Vector3d A2   = MapVector(GetValue(BASE_A2));

    // const Vector3d A3 = Rod_Lam * A03.transpose();
    // const Vector3d A3_1 = Rod_1_Lam * A03.transpose() + Rod_Lam_1 * A03.transpose();
    const Vector3d A3   = MapVector(GetValue(BASE_A3));


    // const double B2 = A2_1.dot(A1);
    // const double B3 = A3_1.dot(A1);

    // const double C12 = A3_1.dot(A2);
    // const double C13 = A2_1.dot(A3);

    // const double Tm = 1 / A11;
    // const double Ts = 1 / A;

    // actual configuration

    const auto phi = ComputeActValue(DISPLACEMENT_ROTATION, 0, shape_functions, GetGeometry());
    const auto phi_1 = ComputeActValue(DISPLACEMENT_ROTATION, 1, shape_functions, GetGeometry());

    const auto x = ComputeActBaseVector(0, shape_functions, GetGeometry());
    const auto a1 = ComputeActBaseVector(1, shape_functions, GetGeometry());
    const auto a1_1 = ComputeActBaseVector(2, shape_functions, GetGeometry());

    const auto a11 = a1.dot(a1);
    const auto a = sqrt(a11);

    const auto t = a1 / a;
    const auto t_1 = a1_1 / a - a1.dot(a1_1) * a1 / pow(a, 3);

    const auto rod = ComputeRod<HyperDual>(t, phi);
    const auto rod_1 = ComputeRod_1<HyperDual>(t, t_1, phi, phi_1);

    const auto lam = ComputeLam<HyperDual>(T, t);
    // const auto lam_1 = ComputeLam_1<HyperDual>(T, T_1, t, t_1);

    const auto rod_lam = rod * lam;
    // const auto rod_1_lam = rod_1 * lam;
    // const auto rod_lam_1 = rod * lam_1;

    // const auto xform = rod_1_lam * Rod_Lam + rod_lam_1 * Rod_Lam + rod_lam * Rod_1_Lam + rod_lam * Rod_Lam_1;

    // const auto y = rod_lam * A2.transpose();
    const Eigen::Matrix<class HyperJet::HyperJet<double>,1,3,1,1,3> a2 = rod_lam * A2.transpose();
    // const auto a2_1 = xform * A02.transpose();

    // const auto z = rod_lam * A3.transpose();
    const Eigen::Matrix<class HyperJet::HyperJet<double>,1,3,1,1,3> a3 = rod_lam * A3.transpose();
    // const auto a3_1 = xform * A03.transpose();

    // const auto b2 = a2_1.dot(a1);
    // const auto b3 = a3_1.dot(a1);

    // const auto c12 = a3_1.dot(a2);
    // const auto c13 = a2_1.dot(a3);

    // Normieren (wobei A2 / A3 eigentlich schon normiert sein sollten!)
    // const auto A22 = A2.dot(A2);
    // const auto A22_length = sqrt(A22);
    // const auto N = A2 / A22_length;

    // const auto A33 = A3.dot(A3);
    // const auto A33_length = sqrt(A33);
    // const auto V = A3 / A33_length;

    const auto d_T = A1.dot(T);
    const auto d_N = A1.dot(A2);
    const auto d_V = A1.dot(A3);

    const auto d_t = a1.dot(T);     // Anteil T in a1 Richtung
    const auto d_n = a1.dot(A2);     // Anteil N in a1 Richtung
    const auto d_v = a1.dot(A3);     // Anteil V in a1 Richtung

    const auto alpha_12 = HyperJet::atan2(a2.dot(A3) , a2.dot(A2));    // Winkel zwischen a2 und A
    const auto alpha_13 = HyperJet::atan2(a3.dot(A2) , a3.dot(A3));    // Winkel zwischen a2 und A
    const auto alpha_2  = HyperJet::atan2(d_n , d_t);                  // Winkel zwischen a1 und A1 um n
    const auto alpha_3  = HyperJet::atan2(d_v , d_t);                  // Winkel zwischen a1 und A1 um v

    const auto u_xyz = x - X;
    const auto delta_a = a1 - A1; 

    // inner energy 
    // auto const condition_type = GetValue(DIRICHLET_CONDITION_TYPE);

    auto const penalty_disp_u   = GetValue(PENALTY_DISPLACEMENT_X);
    auto const penalty_disp_v   = GetValue(PENALTY_DISPLACEMENT_Y);
    auto const penalty_disp_w   = GetValue(PENALTY_DISPLACEMENT_Z);
    auto const penalty_rot    = GetValue(PENALTY_ROTATION);
    auto const penalty_tors   = GetValue(PENALTY_TORSION);


    // std::cout <<"a1 : " << typeid(a1).name() << std::endl;
    // std::cout <<"A1 : " << typeid(A1).name() << std::endl;
    // std::cout <<"a2 : " << typeid(a2).name() << std::endl;
    // std::cout <<"A2 : " << typeid(A2).name() << std::endl;
    // std::cout <<"a3 : " << typeid(a3).name() << std::endl;
    // std::cout <<"A3 : " << typeid(A3).name() << std::endl;
    const auto delta1 = d_t - d_T;
    const auto delta2 = d_n - d_N;
    const auto delta3 = d_v - d_V;

    // Potential Displacement
    // const auto dP_disp = 0.5 * ((delta1).dot(delta1) * pow(penalty_disp_u, 2) 
    //                            +(delta2).dot(delta2) * pow(penalty_disp_v, 2)
    //                            +(delta3).dot(delta3) * pow(penalty_disp_w, 2)) * integration_weight ;

    const auto dP_disp_u = 0.5 * pow((x(0)-X(0)), 2) * penalty_disp_u * integration_weight;                               
    const auto dP_disp_w = 0.5 * pow((x(1)-X(1)), 2) * penalty_disp_w * integration_weight;                               
    const auto dP_disp_v = 0.5 * pow((x(2)-X(2)), 2) * penalty_disp_v * integration_weight;                               

    // const auto dP_disp = 0.5 * ( u_xyz.dot(u_xyz)) * penalty_disp_u * integration_weight ;

    // Potential Rotation
      // By the Gradient 
    const auto dP_grad_bend = 0.5 * (delta_a.dot(delta_a) * penalty_rot) * integration_weight;
      // By the Angle
    const auto dP_alpha_bend = 0.5 * (alpha_2 * alpha_2 + alpha_3 * alpha_3) * penalty_rot * integration_weight; 

    // Potential Torsion
    const auto dP_alpha_tors = 0.5 * (alpha_12 * alpha_12 + alpha_13 * alpha_13) * penalty_tors * integration_weight;
      // const auto dP_alpha_tors =  0.5 * (alpha_12 * alpha_13) * penalty_tors * integration_weight;


    auto const condition_type = 123;

    // Variation Hesse
    if (condition_type == 1)        // Verschiebung 
       MapMatrix(rLeftHandSideMatrix) =    dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h();
    else if (condition_type == 2)    // Torsion
       MapMatrix(rLeftHandSideMatrix) =    dP_alpha_tors.h() ;
    else if (condition_type == 3)   // Biegung  Ableitung
       MapMatrix(rLeftHandSideMatrix) =    dP_alpha_bend.h() ;
    else if (condition_type == 4)   // Biegung Tangens
       MapMatrix(rLeftHandSideMatrix) =    dP_grad_bend.h() ;

    else if(condition_type == 12)   // Verschiebung + Torsion
      MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h()+ dP_alpha_tors.h()) ;
    else if(condition_type == 13)   // Verschiebung + Rotation
      MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h()+ dP_alpha_bend.h()) ;
    else if(condition_type == 23)   // Torsion + Rotation
      MapMatrix(rLeftHandSideMatrix) = ( dP_alpha_tors.h() + dP_alpha_bend.h()) ;
    else if(condition_type == 123)  // Verschiebung + Torsion +  Rotation
      MapMatrix(rLeftHandSideMatrix) = ( dP_disp_u.h() + dP_disp_v.h() + dP_disp_w.h()+ dP_alpha_tors.h() + dP_grad_bend.h()) ;


    // Variation Gradient
    if (condition_type == 1)         // Verschiebung
      MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g()) ;
    else if (condition_type == 2)   // Torsion
      MapVector(rRightHandSideVector) = -( dP_alpha_tors.g()   );
    else if (condition_type == 3)   // Biegung Ableitung
      MapVector(rRightHandSideVector) = -( dP_alpha_bend.g()   );
    else if (condition_type == 4)   // Biegung Tangens
      MapVector(rRightHandSideVector) = -( dP_grad_bend.g() );
    
    else if(condition_type == 12)   // Verschiebung + Torsion
      MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g()+ dP_alpha_tors.g()) ;
    else if(condition_type == 13)   // Verschiebung + Rotation
      MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g()+ dP_alpha_bend.g()) ;
    else if(condition_type == 23)   // Torsion + Rotation
      MapVector(rRightHandSideVector) = -( dP_alpha_tors.g() + dP_alpha_bend.g()) ;
    else if(condition_type == 123)  // Verschiebung + Torsion +  Rotation
      MapVector(rRightHandSideVector) = -( dP_disp_u.g() + dP_disp_v.g() + dP_disp_w.g()+ dP_alpha_tors.g() + dP_grad_bend.g()) ;

    KRATOS_CATCH("")
}

void IgaBeamWeakDirichletCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaBeamWeakDirichletCondition\" #" << Id();
}

} // namespace Kratos