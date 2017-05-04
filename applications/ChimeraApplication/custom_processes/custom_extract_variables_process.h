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

#if !defined(KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED )
#define  KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED


// System includes

#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
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

class CustomExtractVariablesProcess {
public:

	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomExtractVariablesProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomExtractVariablesProcess);

	///@}
	///@name Life Cycle
	///@{

	CustomExtractVariablesProcess(){
		this->pBinLocator = NULL;
		this->pExtractor = NULL;
	}

	/// Destructor.
	virtual ~CustomExtractVariablesProcess() {
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

	template<class TVariableType>
	void ExtractVariable(ModelPart& fromVolumeModelPart, ModelPart& toSurfaceModelPart, Variable<TVariableType>& rVariable){
		this->pExtractor = new CalculateSignedDistanceTo3DConditionSkinProcess(toSurfaceModelPart, fromVolumeModelPart);
		this->pBinLocator = new BinBasedFastPointLocator<3>(fromVolumeModelPart);

		/*
		 * This part of the code below is adapted from "MappingPressureToStructure" function of class CalculateSignedDistanceTo3DSkinProcess
		 */

		{
			//loop over nodes and find the tetra in which it falls, than do interpolation
			array_1d<double, 4 > N;
			const int max_results = 10000;
			BinBasedFastPointLocator<3>::ResultContainerType results(max_results);
			const int n_structure_nodes = toSurfaceModelPart.Nodes().size();

			#pragma omp parallel for firstprivate(results,N)
			//MY NEW LOOP: reset the viisted flaf
			for (int i = 0; i < n_structure_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = toSurfaceModelPart.NodesBegin() + i;
				Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
				p_structure_node->Set(VISITED, false);
			}
			for (int i = 0; i < n_structure_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = toSurfaceModelPart.NodesBegin() + i;
				Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
				BinBasedFastPointLocator<3>::ResultIteratorType result_begin = results.begin();
				Element::Pointer pElement;

				bool is_found = this->pBinLocator->FindPointOnMesh(p_structure_node->Coordinates(), N, pElement, result_begin, max_results);

				if (is_found == true)
				{
					// TODO : For now it only does velocities by components
					// This should be extended to a general variable
					TVariableType nodalVariableValues(0.0);

					Geometry<Node<3> >& geom = pElement->GetGeometry();

					// Do mapping
					ComputeContinuousInterpolation<TVariableType>((*p_structure_node),pElement->GetGeometry(),rVariable,nodalVariableValues);
					//nodalVariableValues is passed by reference

					// Assign Interpolated value to the node
					p_structure_node->FastGetSolutionStepValue(rVariable) = nodalVariableValues;
					p_structure_node->Set(VISITED);
				}
			}
			//AND NOW WE "TREAT" the bad nodes, the ones that belong to the structural faces that by some chance did not cross the fluid elements
			//to such nodes we simply extrapolate the pressure from the neighbors
			int n_bad_nodes=0;
			for (int i = 0; i < n_structure_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = toSurfaceModelPart.NodesBegin() + i;
				Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
				if (p_structure_node->IsNot(VISITED))
					n_bad_nodes++;
			}

			while (n_bad_nodes >= 1.0) {
				int n_bad_nodes_backup = n_bad_nodes;

				for (int i = 0; i < n_structure_nodes; i++) {
					ModelPart::NodesContainerType::iterator iparticle = toSurfaceModelPart.NodesBegin() + i;
					Node < 3 > ::Pointer p_structure_node = *(iparticle.base());

					//here we store the number of neigbor nodes that were given the pressure in the previous loop (i.e. were found)
					if (p_structure_node->IsNot(VISITED)) {
						int n_good_neighbors = 0;
						TVariableType pos_pres(0.0);
						WeakPointerVector< Node < 3 > >& neighours = p_structure_node->GetValue(NEIGHBOUR_NODES);

						for (WeakPointerVector< Node < 3 > >::iterator j = neighours.begin(); j != neighours.end(); j++) {
							if (j->Is(VISITED)) {
								n_good_neighbors++;
								pos_pres += j->FastGetSolutionStepValue(rVariable);
							}
						}
						if (n_good_neighbors != 0) {
							pos_pres /= n_good_neighbors;
							p_structure_node->FastGetSolutionStepValue(rVariable) = pos_pres;
							p_structure_node->Set(VISITED);
							n_bad_nodes--;
						}
					}
				}

				if(n_bad_nodes == n_bad_nodes_backup) break; //WE BREAK THE WHILE HERE, OTHERWISE THE CODE HANGS (it was not able to remove any other node)

			}

		}




		delete pExtractor;
		delete pBinLocator;
	}


	template<class TVariableType>
	void ComputeContinuousInterpolation(const Node<3>& pNode,
            Geometry< Node<3> >& geom,
			Variable<TVariableType>& rVariable,
			TVariableType &nodalVariableValues){

	}


	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	virtual std::string Info() const {
		return "CustomExtractVariablesProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << "CustomExtractVariablesProcess";
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
	CalculateSignedDistanceTo3DConditionSkinProcess* pExtractor;
	//ModelPart &mrBackGroundModelPart;
	//ModelPart &mrPatchSurfaceModelPart;
	BinBasedFastPointLocator<3>* pBinLocator; // Template argument 3 stands for 3D case
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
	CustomExtractVariablesProcess& operator=(CustomExtractVariablesProcess const& rOther);

	/// Copy constructor.
	//CustomExtractVariablesProcess(CustomExtractVariablesProcess const& rOther);

	///@}

}; // Class CustomExtractVariablesProcess

}  // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
