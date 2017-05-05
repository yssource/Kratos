// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================
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
#include "custom_processes/apply_multi_point_constraints_process_chimera.h"
// Application includes
#include "custom_utilities/multipoint_constraint_data.hpp"

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

	CustomApplyMpcConstraint2dProcess(ModelPart& surfaceModelPart){		
		this->pBinLocator = BinBasedFastPointLocator<2>::Pointer( new BinBasedFastPointLocator<2>(surfaceModelPart) );
		this->pMpcProcess = ApplyMultipointConstraintsProcessChimera::Pointer( new ApplyMultipointConstraintsProcessChimera(surfaceModelPart) );		
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

	
	void ApplyMpcConstraint2d(ModelPart& boundaryModelPart){

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
					

					Geometry<Node<3> >& geom = pElement->GetGeometry();
					//std::cout<<"Element Id: "<< pElement->Id();

					// TO DO Apply MPC constraints. Here have to use the mpc functions to apply the constraint
					
					//std::cout<<p_boundary_node->Id()<<" \t"<<geom[0].Id()<<" \t"<<geom[1].Id()<<" \t"<<geom[2].Id()<<" \t"<<N[0]<<" \t "<<N[1]<<" \t"<< N[2]<<std::endl;										
					
					for(int i = 0; i < geom.size(); i++){
					
					pMpcProcess->AddMasterSlaveRelation( geom[i],VELOCITY_X,*p_boundary_node,VELOCITY_X,N[i], 0 );
					pMpcProcess->AddMasterSlaveRelation( geom[i],VELOCITY_Y,*p_boundary_node,VELOCITY_Y,N[i], 0 );
				
					
					}	//Nodes loop inside the found element ends here
										
					
					p_boundary_node->Set(VISITED);
					
				}
			}
			
	
			std::cout<<"mpcProcess ends"<<std::endl;

		}
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
	BinBasedFastPointLocator<2>::Pointer pBinLocator; // Template argument 3 stands for 3D case
	ApplyMultipointConstraintsProcessChimera::Pointer pMpcProcess;
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
