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
#include <iostream>
#include <iomanip>
#include "includes/define.h"
#include "includes/variables.h"
#include <chrono>
// External includes

// Project includes
#include "iga_base_condition.h"
#include "iga_beam_load_condition.h"
#include "iga_application_variables.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	IgaBeamLoadCondition::IgaBeamLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: IgaBaseCondition<4>(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	IgaBeamLoadCondition::IgaBeamLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: IgaBaseCondition<4>(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer IgaBeamLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return IgaBaseCondition<4>::Pointer(new IgaBeamLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	IgaBeamLoadCondition::~IgaBeamLoadCondition()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void IgaBeamLoadCondition::EquationIdVector(
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

	//************************************************************************************
	//************************************************************************************
	void IgaBeamLoadCondition::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
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

	void IgaBeamLoadCondition::CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool ComputeLeftHandSide,
		const bool ComputeRightHandSide)
	{
		KRATOS_TRY
		using namespace BeamUtilities;
		using std::sqrt;
		using std::pow;

		using Vector3d = BeamUtilities::Vector<3>;
		using Matrix3d = BeamUtilities::Matrix<3, 3>;

		// get integration weight

		// const double& integration_weight = GetValue(INTEGRATION_WEIGHT);

		// get shape functions

		BeamUtilities::Matrix<3, Dynamic> shape_functions(3, GetValue(SHAPE_FUNCTION_VALUES).size());
		shape_functions.row(0) = MapVector(this->GetValue(SHAPE_FUNCTION_VALUES));
		shape_functions.row(1) = MapVector(this->GetValue(SHAPE_FUNCTION_LOCAL_DER_1));
		shape_functions.row(2) = MapVector(this->GetValue(SHAPE_FUNCTION_LOCAL_DER_2));

		const Vector3d A1   = MapVector(this->GetValue(BASE_A1));
		const Vector3d A1_1 = MapVector(this->GetValue(BASE_A1_1));

		const double A11 = A1.dot(A1);
		const double A = sqrt(A11);

		const Vector3d T = A1 / A;
		const Vector3d T_1 = A1_1 / A - A1.dot(A1_1) * A1 / pow(A, 3);

		const Vector3d A2   = MapVector(GetValue(BASE_A2));
		const Vector3d A3   = MapVector(GetValue(BASE_A3));

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
		const auto lam_1 = ComputeLam_1<HyperDual>(T, T_1, t, t_1);

		const auto rod_lam = rod * lam;

		const auto a2 = rod_lam * A2.transpose();
		const auto a3 = rod_lam * A3.transpose();


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

		array_1d<double, 3 > PointLoad = ZeroVector(3);

		PointLoad = GetValue(LOAD_VECTOR_MOMENT); 

		auto const dP_x = 0.5 * (alpha_12 + alpha_13) * PointLoad[0];
		auto const dP_z = alpha_2 * PointLoad[2]; 
		auto const dP_y = alpha_3 * PointLoad[1];
		
		MapMatrix(rLeftHandSideMatrix) =  dP_y.h() + dP_z.h() ;
		MapVector(rRightHandSideVector) = -(dP_y.g() + dP_z.g());

		KRATOS_CATCH( "" )
    }


} // Namespace Kratos
