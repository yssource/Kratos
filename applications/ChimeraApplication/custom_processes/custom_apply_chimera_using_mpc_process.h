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

#if !defined(KRATOS_CUSTOM_APPLY_CHIMERA_USING_CHIMERA_H_INCLUDED)
#define KRATOS_CUSTOM_APPLY_CHIMERA_USING_CHIMERA_H_INCLUDED

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
#include "includes/variables.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

// Application includes
#include "custom_utilities/multipoint_constraint_data.hpp"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_hole_cutting_process.h"
#include "custom_processes/apply_multi_point_constraints_process.h"

namespace Kratos
{

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

template <unsigned int TDim>

class CustomApplyChimeraUsingMpcProcess
{
  public:
	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomExtractVariablesProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomApplyChimeraUsingMpcProcess);

	typedef typename BinBasedFastPointLocator<TDim>::Pointer BinBasedPointLocatorPointerType;
	///@}
	///@name Life Cycle
	///@{
	CustomApplyChimeraUsingMpcProcess(ModelPart &AllModelPart, ModelPart &BackgroundModelPart, ModelPart &PatchModelPart, double distance = 1e-12) : mrAllModelPart(AllModelPart), mrBackgroundModelPart(BackgroundModelPart), mrPatchModelPart(PatchModelPart), overlap_distance(distance)

	{
		this->pBinLocatorForBackground = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(mrBackgroundModelPart));
		this->pBinLocatorForPatch = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(mrPatchModelPart));
		this->pMpcProcess = ApplyMultipointConstraintsProcess::Pointer(new ApplyMultipointConstraintsProcess(mrAllModelPart));
		this->pHoleCuttingProcess = CustomHoleCuttingProcess::Pointer(new CustomHoleCuttingProcess());
		this->pCalculateDistanceProcess = typename CustomCalculateSignedDistanceProcess<TDim>::Pointer(new CustomCalculateSignedDistanceProcess<TDim>());
	}

	/// Destructor.
	virtual ~CustomApplyChimeraUsingMpcProcess()
	{
	}

	///@}
	///@name Operators
	///@{

	void operator()()
	{
		Execute();
	}

	///@}
	///@name Operations
	///@{

	virtual void Execute()
	{
	}

	virtual void Clear()
	{
	}

	void ApplyMpcConstraint(ModelPart &mrBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator, unsigned int type = 1)
	{

		{
			//loop over nodes and find the triangle in which it falls, than do interpolation
			array_1d<double, TDim + 1> N;
			const int max_results = 10000;
			typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
			const int n_boundary_nodes = mrBoundaryModelPart.Nodes().size();

#pragma omp parallel for firstprivate(results, N)
			//MY NEW LOOP: reset the visited flag
			for (int i = 0; i < n_boundary_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = mrBoundaryModelPart.NodesBegin() + i;
				Node<3>::Pointer p_boundary_node = *(iparticle.base());
				p_boundary_node->Set(VISITED, false);
			}

			for (int i = 0; i < n_boundary_nodes; i++)

			{
				ModelPart::NodesContainerType::iterator iparticle = mrBoundaryModelPart.NodesBegin() + i;
				Node<3>::Pointer p_boundary_node = *(iparticle.base());

				typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

				Element::Pointer pElement;

				bool is_found = false;
				is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

				if (is_found == true)
				{
					// TODO : For now it only does velocities by components
					// This should be extended to a general variable

					//Check if some of the host elements are made inactive

					if ((pElement)->IsDefined(ACTIVE))
					{

						if (!(pElement->Is(ACTIVE)))
						{
							std::cout << "Warning : One of the hole element is used for MPC constraint" << std::endl;
							pElement->Set(ACTIVE);
							std::cout << "Setting the element: " << pElement->Id() << " to active" << std::endl;
						}
					}

					

					Geometry<Node<3>> &geom = pElement->GetGeometry();

					{

						for (int i = 0; i < geom.size(); i++)
						{

							pMpcProcess->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_X, *p_boundary_node, VELOCITY_X, N[i]);
							pMpcProcess->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_Y, *p_boundary_node, VELOCITY_Y, N[i]);

							if (TDim == 3)
								pMpcProcess->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_Z,*p_boundary_node, VELOCITY_Z, N[i]);
							//pMpcProcess->AddMasterSlaveRelationVariables( geom[i],PRESSURE,*p_boundary_node,PRESSURE,N[i], 0 );
						}
					}
				}
			}

			if (type == 0)

			{
				ModelPart::NodesContainerType::iterator iparticle = mrBoundaryModelPart.NodesBegin();
				Node<3>::Pointer p_boundary_node = *(iparticle.base());

				typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

				Element::Pointer pElement;

				bool is_found = false;
				is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

				if (is_found == true)
				{

					Geometry<Node<3>> &geom = pElement->GetGeometry();

					std::cout << "Presssure of " << p_boundary_node->Id() << "coupled to " << pElement->Id() << std::endl;

					for (int i = 0; i < geom.size(); i++)
					{

						pMpcProcess->AddMasterSlaveRelationWithNodesAndVariable(geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);
					}
				}

				else
				{
					std::cout << "Cannot find host element for Pressure coupling" << std::endl;
					std::exit(-1);
				}
			}

			std::cout << "mpcProcess ends" << std::endl;
		}
	}

	//Apply Chimera with or without overlap
	void ApplyChimeraUsingMpc(ModelPart &mrPatchBoundaryModelPart)
	{
		const double epsilon = 1e-12;
		if (overlap_distance < epsilon)
		{
			std::cout << "Overlap distance should be a positive number" << std::endl;
			std::exit(-1);
		}
		if (overlap_distance > epsilon)

		{

			ModelPart HoleModelPart = ModelPart("HoleModelpart");
			ModelPart HoleBoundaryModelPart = ModelPart("HoleBoundaryModelPart");

			//this->pCalculateDistanceProcess->ExtractDistance(mrPatchModelPart, mrBackgroundModelPart, mrPatchBoundaryModelPart);
			this->pCalculateDistanceProcess->CalculateSignedDistance(mrBackgroundModelPart, mrPatchBoundaryModelPart);

			this->pHoleCuttingProcess->CreateHoleAfterDistance(mrBackgroundModelPart, HoleModelPart, HoleBoundaryModelPart, overlap_distance);

			ApplyMpcConstraint(mrPatchBoundaryModelPart, pBinLocatorForBackground, 0);
			std::cout << "Patch boundary coupled with background" << std::endl;
			ApplyMpcConstraint(HoleBoundaryModelPart, pBinLocatorForPatch, 1);
			std::cout << "HoleBoundary  coupled with patch" << std::endl;
		}

		else
		{

			std::cout << "Applying Strong MPC " << std::endl;

			ModelPart HoleModelPart = ModelPart("HoleModelPart");
			ModelPart HoleBoundaryModelPart = ModelPart("HoleBoundaryModelPart");

			//this->pCalculateDistanceProcess->ExtractDistance(mrPatchModelPart, mrBackgroundModelPart, mrPatchBoundaryModelPart);
			this->pCalculateDistanceProcess->CalculateSignedDistance(mrBackgroundModelPart, mrPatchBoundaryModelPart);

			this->pHoleCuttingProcess->CreateHoleAfterDistance(mrBackgroundModelPart, HoleModelPart, HoleBoundaryModelPart, epsilon);

			ApplyMpcConstraint(mrPatchBoundaryModelPart, pBinLocatorForBackground, 1);
		}
	}

	void SetOverlapDistance(double distance)
	{

		this->overlap_distance = distance;
	}

	virtual std::string Info() const
	{
		return "CustomApplyChimeraUsingMpcProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "CustomApplyChimeraUsingMpcProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
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
	BinBasedPointLocatorPointerType pBinLocatorForBackground; // Template argument 3 stands for 3D case
	BinBasedPointLocatorPointerType pBinLocatorForPatch;
	ApplyMultipointConstraintsProcess::Pointer pMpcProcess;
	CustomHoleCuttingProcess::Pointer pHoleCuttingProcess;
	typename CustomCalculateSignedDistanceProcess<TDim>::Pointer pCalculateDistanceProcess;
	ModelPart &mrAllModelPart;
	ModelPart &mrBackgroundModelPart;
	ModelPart &mrPatchModelPart;
	double overlap_distance;

	// epsilon
	//static const double epsilon;

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
	CustomApplyChimeraUsingMpcProcess &operator=(CustomApplyChimeraUsingMpcProcess const &rOther);

	/// Copy constructor.
	//CustomExtractVariablesProcess(CustomExtractVariablesProcess const& rOther);

	///@}

}; // Class CustomExtractVariablesProcess

} // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
