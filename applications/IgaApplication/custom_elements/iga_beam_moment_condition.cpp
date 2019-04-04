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
#include "iga_beam_moment_condition.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaBeamMomentCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaBeamMomentCondition>(NewId, geometry,
        pProperties);
}

void IgaBeamMomentCondition::GetDofList(
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

void IgaBeamMomentCondition::EquationIdVector(
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

void IgaBeamMomentCondition::Initialize()
{
    
}

void IgaBeamMomentCondition::CalculateAll(
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
    const double prestress = properties[PRESTRESS_CAUCHY];
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

    // const Vector3d A1 = ComputeRefBaseVector(1, shape_functions, GetGeometry());
    // const Vector3d A1_1 = ComputeRefBaseVector(2, shape_functions, GetGeometry());
    const Vector3d A1   = MapVector(GetValue(BASE_A1));
    const Vector3d A1_1 = MapVector(GetValue(BASE_A1_1));

    const double A11 = A1.dot(A1);
    const double A = sqrt(A11);

    const Vector3d T = A1 / A;
    const Vector3d T_1 = A1_1 / A - A1.dot(A1_1) * A1 / pow(A, 3);

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
    const Vector3d A2_1 = MapVector(GetValue(BASE_A2_1));

    // const Vector3d A3 = Rod_Lam * A03.transpose();
    // const Vector3d A3_1 = Rod_1_Lam * A03.transpose() + Rod_Lam_1 * A03.transpose();
    // const Vector3d A3   = T.cross(A2); 
    // const Vector3d A3_1 = T.cross(A2_1) + T_1.cross(A2); 
    const Vector3d A3   = MapVector(GetValue(BASE_A3));
    const Vector3d A3_1 = MapVector(GetValue(BASE_A3_1));

    const double B2 = A2_1.dot(A1);
    const double B3 = A3_1.dot(A1);

    const double C12 = A3_1.dot(A2);
    const double C13 = A2_1.dot(A3);

    const double Tm = 1 / A11;
    const double Ts = 1 / A;

    // actual configuration

    const auto phi = ComputeActValue(DISPLACEMENT_ROTATION, 0, shape_functions, GetGeometry());
    const auto phi_1 = ComputeActValue(DISPLACEMENT_ROTATION, 1, shape_functions, GetGeometry());    

    const auto a1 = ComputeActBaseVector(1, shape_functions, GetGeometry());
    const auto a1_1 = ComputeActBaseVector(2, shape_functions, GetGeometry());

    const auto a11 = a1.dot(a1);
    const auto a = sqrt(a11);

    const auto t = a1 / a;
    const auto t_1 = a1_1 / a - a1.dot(a1_1) * a1 / pow(a, 3);

    const auto rod = ComputeRod<HyperDual>(t, phi);
    const auto rod_1 = ComputeRod_1<HyperDual>(t, t_1, phi, phi_1);

    const auto lam = ComputeLam<HyperDual>(T, t);
    const auto lam_1 = ComputeLam_1<HyperDual>(T, T_1, t, t_1);

    const auto rod_lam = rod * lam;
    const auto rod_1_lam = rod_1 * lam;
    const auto rod_lam_1 = rod * lam_1;

    // const auto xform = rod_1_lam * Rod_Lam + rod_lam_1 * Rod_Lam + rod_lam * Rod_1_Lam + rod_lam * Rod_Lam_1;

    const auto a2 = rod_lam * A2.transpose();
    // const auto a2_1 = xform * A02.transpose();
    const auto a2_1 = rod_lam * A2_1.transpose() + rod_lam_1 * A2.transpose() + rod_1_lam * A2.transpose() ;

    const auto a3 = rod_lam * A3.transpose();
    // const auto a3_1 = xform * A03.transpose();
    const auto a3_1 = rod_lam * A3_1.transpose() + rod_lam_1 * A3.transpose() + rod_1_lam * A3.transpose();

    const auto b2 = a2_1.dot(a1);
    const auto b3 = a3_1.dot(a1);

    const auto c12 = a3_1.dot(a2);
    const auto c13 = a2_1.dot(a3);
    
    const auto _load_vec = MapVector(GetValue(LOAD_VECTOR_MOMENT));

    const auto fac_t = a1.dot(_load_vec); 
    const auto fac_n = a3.dot(_load_vec);
    const auto fac_v = a2.dot(_load_vec);

    // inner energy

    const auto eps11 = Tm * (a11 - A11) / 2;
    const auto kap2  = Tm * (b2  - B2 );
    const auto kap3  = Tm * (b3  - B3 );
    const auto kap12 = Ts * (c12 - C12);
    const auto kap13 = Ts * (c13 - C13);

    const auto dP = 0.5 * (pow(kap2 , 2)* ei3 * fac_v +
                           pow(kap3 , 2)* ei2 * fac_n) * A * integration_weight;

    // // normalize base vectors
    // Vector3 t = r1 / norm_2(r1);
    // n /= norm_2(n);
    // v /= norm_2(v);

    // double fac_t;
    // double fac_n;
    // double fac_v;

    // fac_t = inner_prod(t, _load_vec);        //TODO check if right length
    // fac_n = inner_prod(n, _load_vec);
    // fac_v = inner_prod(v, _load_vec);
    // // LOG("_load_vec " << _load_vec);
    // // LOG("fac_t " << fac_t);
    // // LOG("fac_n " << fac_n);
    // // LOG("fac_v " << fac_v);


    // for(size_t k = 0; k != NumberOfDofs(); k++)
    // {
    //     _local_load_vec[k] = curv_var_t[k] * fac_t + curv_var_n[k] * fac_n + curv_var_v[k] * fac_v;
    //     // _local_load_vec[k] = 0;
    // }



    MapMatrix(rLeftHandSideMatrix) = dP.h();
    MapVector(rRightHandSideVector) = -dP.g();

    KRATOS_CATCH("");
}

void IgaBeamMomentCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaBeamMomentCondition\" #" << Id();
}

} // namespace Kratos





// // LOADS________________________________________________________________________________________________

// void IgaBeamMomentCondition::CalculateLoadMoment(
//     VectorType& _local_load_vec,
//         Vector3 _load_vec)
// {
// KRATOS_TRY;

//     // resize loadvector
//     _local_load_vec .resize(NumberOfDofs());
//     _local_load_vec.clear();

//     // Get Shape Functions
//     // Vector& shape_function_values = GetValue(SHAPE_FUNCTION_VALUES);
//     // Matrix& shape_derivatives     = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
//     //Get Tangent Vector
//     Vector3 T0_vec = GetValue(T0);
//     // Normalize t0
//     double t0_L = norm_2(T0_vec);
//     T0_vec = T0_vec/t0_L;
//     // get Rotations
//     double Phi      = GetValue(PHI);
//     double Phi_der  = GetValue(PHI_DER_1);
//     double phi      = 0;
//     double phi_der  = 0;
    
//     // Declarations
//     Vector3 R1;
//     Vector3 R2;
//     double A;
//     double B;
//     Vector3 r1;
//     Vector3 r2;
//     double a;
//     double b;
//     R1.clear();
//     R2.clear();
//     r1.clear();
//     r2.clear();
    
//     double B_n;
//     double B_v;
//     double C_12;
//     double C_13;
//     double b_n;
//     double b_v;
//     double c_12;
//     double c_13;

//     Vector3 N;     // Principal Axis 1 of Cross Section in Undeformed Config.
//     Vector3 V;   // Prinzipal Axis 2 of Cross Section in Undeformed Config.
//     Vector3 n;   // Principal Axis 1 Of Cross Section in Deformed Config.
//     Vector3 v;  // Principal Axis 2 of Cross Section in Deformed Config.
//     Vector3 N0;    // Principal Axis 1 of Cross Section in Reference Config.
//     Vector3 V0;    // Principal Axis 2 of Cross Section in Reference Config.

//     VectorType curv_var_t;
//     VectorType curv_var_n;
//     VectorType curv_var_v;

//    // Compute the Vectors R1 R2 and the length A and B in undeformed and deformed state
//     ComputeGeometryReference(R1, R2, A, B);
//     ComputeGeometryActual(r1, r2, a, b);
//     ComputeCrossSectionGeometryReference(R1, R2, N, V, T0_vec, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der);
//     ComputeCrossSectionGeometryActual(R1, R2, r1, r2, N0, V0, n, v, b_n, b_v, c_12, c_13, Phi, Phi_der, phi, phi_der);
//     ComputePhiReferenceProperty(phi, phi_der);
//     ComputeRotationalDof( curv_var_t, curv_var_n, curv_var_v, R1, r1, N0, V0, N, V, n, v, Phi, phi );
//     // LOG("R1" << R1);
//     // LOG("r1" << r1);
//     // LOG("r2" << r2);


//     // normalize base vectors
//     Vector3 t = r1 / norm_2(r1);
//     n /= norm_2(n);
//     v /= norm_2(v);

//     double fac_t;
//     double fac_n;
//     double fac_v;

//     fac_t = inner_prod(t, _load_vec);        //TODO check if right length
//     fac_n = inner_prod(n, _load_vec);
//     fac_v = inner_prod(v, _load_vec);
//     // LOG("_load_vec " << _load_vec);
//     // LOG("fac_t " << fac_t);
//     // LOG("fac_n " << fac_n);
//     // LOG("fac_v " << fac_v);


//     for(size_t k = 0; k != NumberOfDofs(); k++)
//     {
//         _local_load_vec[k] = curv_var_t[k] * fac_t + curv_var_n[k] * fac_n + curv_var_v[k] * fac_v;
//         // _local_load_vec[k] = 0;
//     }

//     // LOG("curv_var_t : " << curv_var_t); 
//     // LOG("curv_var_n : " << curv_var_n); 
//     // LOG("curv_var_v : " << curv_var_v); 
//     // LOG("Load Vector: " << _local_load_vec); 

// // // get properties
// //     const auto& properties  = GetProperties();
// //     const double moment  = properties[LOAD_MOMENT];        //TODO: könnte eigentlich auch direkt in der Funktion abgerufen werden

// //     _local_load_vec[NumberOfDofs() -2] = moment;


// KRATOS_CATCH(""); 
// }







// void IgaBeamMomentCondition::PrintInfo(std::ostream& rOStream) const
// {
//     rOStream << "\"IgaBeamMomentCondition\" #" << Id();
// }

// } // namespace Kratos