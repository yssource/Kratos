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
#include <chrono>
// External includes

// Project includes
#include "iga_beam_ad_element.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaBeamADElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaBeamADElement>(NewId, geometry,
        pProperties);
}

void IgaBeamADElement::GetDofList(
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

void IgaBeamADElement::EquationIdVector(
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

void IgaBeamADElement::Initialize()
{
    
}

void IgaBeamADElement::CalculateAll(
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


    std::ofstream write;
    write.open("OutputAD.txt", std::ofstream::app);
    auto start = std::chrono::steady_clock::now();

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

    // inner energy

    const auto eps11 = Tm * (a11 - A11) / 2;
    const auto kap2  = Tm * (b2  - B2 );
    const auto kap3  = Tm * (b3  - B3 );
    const auto kap12 = Ts * (c12 - C12);
    const auto kap13 = Ts * (c13 - C13);

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    write << time.count() <<  "\n";
    write.close();

    const auto dP = 0.5 * (pow(eps11, 2) * ea +
                           pow(kap2 , 2) * ei3 +
                           pow(kap3 , 2) * ei2 +
                           pow(kap12, 2) * gi1 * 0.5 +
                           pow(kap13, 2) * gi1 * 0.5) * A * integration_weight;


    MapMatrix(rLeftHandSideMatrix) = dP.h();
    MapVector(rRightHandSideVector) = -dP.g();

    KRATOS_CATCH("")
}

void IgaBeamADElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaBeamADElement\" #" << Id();
}

} // namespace Kratos