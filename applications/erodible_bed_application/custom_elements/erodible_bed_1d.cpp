//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/erodible_bed_1d.h"
#include "erodible_bed_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	ErodibleBed1D::ErodibleBed1D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ErodibleBed1D::ErodibleBed1D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer ErodibleBed1D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ErodibleBed1D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ErodibleBed1D::~ErodibleBed1D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ErodibleBed1D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		boost::numeric::ublas::bounded_matrix<double,2,1> msDN_DX;  //gradients matrix 
		boost::numeric::ublas::bounded_matrix<double,1,1> msD;  //conductivity matrix 
		msD = ZeroMatrix(1,1); //initializing the matrix as zero
  		array_1d<double,2> msN; //dimension = number of nodes . Position of the gauss point 
		array_1d<double,2> ms_temp; //dimension = number of nodes . 
		array_1d<double,2> old_temp; //dimension = number of nodes . 

		const double delta_t = rCurrentProcessInfo[DELTA_TIME];

		const unsigned int number_of_points = GetGeometry().size();
		
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			old_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE,1);

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); //resizing the system in case it does not have the right size 
		rLeftHandSideMatrix=ZeroMatrix(number_of_points,number_of_points);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);
		rRightHandSideVector = ZeroVector(number_of_points);

		//getting data for the given geometry
		const double lenght =GetGeometry().Length();
		msN[0]=0.5;
		msN[1]=0.5;
		msDN_DX(0,0) = 1.0/(GetGeometry()[0].X()-GetGeometry()[1].X());
		msDN_DX(1,0) = 1.0/(GetGeometry()[1].X()-GetGeometry()[0].X());

		double conductivity = 1.0; // GetProperties()[CONDUCTIVITY];
		msD(0,0)=conductivity;
		
		//normal Laplacian matrix
		boost::numeric::ublas::bounded_matrix<double,2,2> Laplacian_matrix = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));  //gradients matrix 
		Laplacian_matrix*=lenght;
		//adding the laplacian to the LHS and rhs
		const double theta=0.5;
		rLeftHandSideMatrix += theta*Laplacian_matrix;  // Bt D B
		noalias(rRightHandSideVector) -= (1.0-theta)*prod(Laplacian_matrix,old_temp);

		//normal mass matrix
		double heat_capacity = 1.0;
		boost::numeric::ublas::bounded_matrix<double,2,2> Mass_matrix=ZeroMatrix(number_of_points,number_of_points);
		Mass_matrix(0,0) = heat_capacity * 0.5;
		Mass_matrix(1,1) = heat_capacity * 0.5;
		Mass_matrix*=lenght;
		//adding the mass matrix to the rhs
		rLeftHandSideMatrix += (1.0/delta_t)*Mass_matrix;
		noalias(rRightHandSideVector) += (1.0/delta_t)*prod(Mass_matrix,old_temp);
		

		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void ErodibleBed1D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void ErodibleBed1D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void ErodibleBed1D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void ErodibleBed1D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}



} // Namespace Kratos
