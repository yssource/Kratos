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
#include "iga_beam_load_condition.h"
#include "iga_application_variables.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	IgaBeamLoadCondition::IgaBeamLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	IgaBeamLoadCondition::IgaBeamLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer IgaBeamLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new IgaBeamLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	IgaBeamLoadCondition::~IgaBeamLoadCondition()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void IgaBeamLoadCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rRightHandSideVector.size() != 1)
			rRightHandSideVector.resize(1,false);

		// double load = GetGeometry()[0].GetSolutionStepValue(POINT_HEAT_SOURCE);
		// rRightHandSideVector[0] = load;
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void IgaBeamLoadCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if(rLeftHandSideMatrix.size1() != 1)
			rLeftHandSideMatrix.resize(1,1,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(1,1);
		if(rRightHandSideVector.size() != 1)
			rRightHandSideVector.resize(1,false);
		// double load = GetGeometry()[0].GetSolutionStepValue(POINT_HEAT_SOURCE);
		// rRightHandSideVector[0] = load;
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void IgaBeamLoadCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(TEMPERATURE).EquationId());			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void IgaBeamLoadCondition::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        // KRATOS_TRY;

        // rElementalDofList.resize(NumberOfDofs());

        // for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        //     SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        //     SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        //     SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
        //     SetElementDof(rElementalDofList, i, 3, DISPLACEMENT_ROTATION);
        // }

        // KRATOS_CATCH("")
    }



} // Namespace Kratos
