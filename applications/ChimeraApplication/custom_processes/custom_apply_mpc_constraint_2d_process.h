//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, , KratosAppGenerator
//

#if !defined(KRATOS_CUSTOM_APPLY_MPC_CONSTRAINT_2D_H_INCLUDED )
#define  KRATOS_CUSTOM_APPLY_MPC_CONSTRAINT_2D_H_INCLUDED


// System includes

#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"

#include "utilities/binbased_fast_point_locator.h"

// Project includes

#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.

class CustomApplyMpcConstraint2dProcess {
public:

	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomExtractVariablesProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomApplyMpcConstraint2dProcess);

	///@}
	///@name Life Cycle
	///@{

	CustomApplyMpcConstraint2dProcess(){
		this->pBinLocator = NULL;
		
	}

	/// Destructor.
	virtual ~CustomApplyMpcConstraint2dProcess() {
	}

	///@}
	///@name Operators
	///@{

	void operator()() {
		Execute();
	}

	///@}
	///@name Operations
	///@{

	virtual void Execute() {
	}

	virtual void Clear() {
	}

	
	void ApplyMpcConstraint2d(ModelPart& surfaceModelPart, ModelPart& boundaryModelPart){
		
		this->pBinLocator = new BinBasedFastPointLocator<2>(surfaceModelPart);

		/*
		 * This part of the code below is adapted from "MappingPressureToStructure" function of class CalculateSignedDistanceTo3DSkinProcess
		 */

		{
			//loop over nodes and find the triangle in which it falls, than do interpolation
			array_1d<double, 3 > N;
			const int max_results = 10000;
			BinBasedFastPointLocator<2>::ResultContainerType results(max_results);
			const int n_boundary_nodes = boundaryModelPart.Nodes().size();
			

			#pragma omp parallel for firstprivate(results,N)
			//MY NEW LOOP: reset the visited flag
			for (int i = 0; i < n_boundary_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = boundaryModelPart.NodesBegin() + i;
				Node < 3 > ::Pointer p_boundary_node = *(iparticle.base());
				p_boundary_node->Set(VISITED, false);
			}
			for (int i = 0; i < n_boundary_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = boundaryModelPart.NodesBegin() + i;
				Node < 3 > ::Pointer p_boundary_node = *(iparticle.base());
				
				BinBasedFastPointLocator<2>::ResultIteratorType result_begin = results.begin();
				
				Element::Pointer pElement;
				
				bool is_found = false;
				is_found = this->pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement,result_begin, max_results);
				
				if (is_found == true)
				{
					// TODO : For now it only does velocities by components
					// This should be extended to a general variable
					

					//Geometry<Node<3> >& geom = pElement->GetGeometry();

					// TO DO Apply MPC constraints. Here have to use the mpc functions to apply the constraint
					//std::cout<<"Node Id: "<<p_boundary_node->Id()<<std::endl;
					
					

					//std::cout<<"Element Id: "<< pElement->Id();
					
										
					//nodalVariableValues is passed by reference
					
					// Assign Interpolated value to the node
					
					p_boundary_node->Set(VISITED);
				}
			}
			//AND NOW WE "TREAT" the bad nodes, the ones that belong to the structural faces that by some chance did not cross the fluid elements
			//to such nodes we simply extrapolate the pressure from the neighbors
	

			

		}




		
		delete pBinLocator;
	}


	virtual std::string Info() const {
		return "CustomApplyMpcConstraint2dProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << "CustomApplyMpcConstraint2dProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {
	}

	///@}
	///@name Friends
	///@{

	///@}

protected:

	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:

	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{
	
	//ModelPart &mrBackGroundModelPart;
	//ModelPart &mrPatchSurfaceModelPart;
	BinBasedFastPointLocator<2>* pBinLocator; // Template argument 3 stands for 3D case
	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	CustomApplyMpcConstraint2dProcess& operator=(CustomApplyMpcConstraint2dProcess const& rOther);

	/// Copy constructor.
	//CustomExtractVariablesProcess(CustomExtractVariablesProcess const& rOther);

	///@}

}; // Class CustomExtractVariablesProcess

}  // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
