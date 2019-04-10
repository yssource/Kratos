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

		// const unsigned int NumberOfNodes = GetGeometry().size();
		// const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

		// // Resizing as needed the LHS
		// const unsigned int MatSize = NumberOfNodes * Dimension;

		// if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
		// {
		// 	if ( rLeftHandSideMatrix.size1() != MatSize )
		// 	{
		// 		rLeftHandSideMatrix.resize( MatSize, MatSize, false );
		// 	}

		// 	noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
		// }

		// //resizing as needed the RHS
		// if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
		// {
		// 	if ( rRightHandSideVector.size( ) != MatSize )
		// 	{
		// 		rRightHandSideVector.resize( MatSize, false );
		// 	}

		// 	noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
		// }

		// // Vector with a loading applied to the condition
		// array_1d<double, 3 > PointLoad = ZeroVector(3);
		// if( this->Has( POINT_LOAD ) )
		// {
		// 	noalias(PointLoad) = this->GetValue( POINT_LOAD );
		// }

		// for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
		// {
		// 	const unsigned int base = ii*Dimension;

		// 	if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
		// 	{
		// 		noalias(PointLoad) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
		// 	}

		// 	for(unsigned int k = 0; k < Dimension; ++k)
		// 	{
		// 		rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * PointLoad[k];
		// 	}
		// }

		KRATOS_CATCH( "" )
    }


} // Namespace Kratos
