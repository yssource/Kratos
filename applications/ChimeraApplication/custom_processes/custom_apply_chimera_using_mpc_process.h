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
#include "utilities/math_utils.h"
# include "includes/kratos_parameters.h"

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
#include "custom_utilities/vtk_output.hpp"

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
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	typedef std::pair<unsigned int, unsigned int> SlavePairType;
	typedef Kratos::MpcData::MasterDofWeightMapType MasterDofWeightMapType;
	typedef ProcessInfo ProcessInfoType;
	typedef MpcData::Pointer MpcDataPointerType;
	typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
	typedef Dof<double> DofType;
	typedef std::vector<DofType> DofVectorType;
	typedef MpcData::VariableComponentType VariableComponentType;

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

	void ApplyMpcConstraint(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator, unsigned int type = 1)
	{

		{
			//loop over nodes and find the triangle in which it falls, than do interpolation
			array_1d<double, TDim + 1> N;
			const int max_results = 10000;
			typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
			const int n_boundary_nodes = rBoundaryModelPart.Nodes().size();

#pragma omp parallel for firstprivate(results, N)
			//MY NEW LOOP: reset the visited flag
			for (int i = 0; i < n_boundary_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin() + i;
				Node<3>::Pointer p_boundary_node = *(iparticle.base());
				p_boundary_node->Set(VISITED, false);
			}

			for (int i = 0; i < n_boundary_nodes; i++)

			{
				ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin() + i;
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

					/*if ((pElement)->IsDefined(ACTIVE))
					{

						if (!(pElement->Is(ACTIVE)))
						{
							std::cout << "Warning : One of the hole element is used for MPC constraint" << std::endl;
							pElement->Set(ACTIVE);
							std::cout << "Setting the element: " << pElement->Id() << " to active" << std::endl;
						}
					}*/

					Geometry<Node<3>> &geom = pElement->GetGeometry();

					{

						for (int i = 0; i < geom.size(); i++)
						{

							pMpcProcess->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_X, *p_boundary_node, VELOCITY_X, N[i]);
							pMpcProcess->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_Y, *p_boundary_node, VELOCITY_Y, N[i]);

							if (TDim == 3)
								pMpcProcess->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_Z, *p_boundary_node, VELOCITY_Z, N[i]);
							pMpcProcess->AddMasterSlaveRelationWithNodesAndVariable(geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);
						}
					}
				}
			}

			if (type == 0)

			{

				//double distance;
				//const int n_patch_nodes = mrPatchModelPart.Nodes().size();
				//bool IsCoupled = false;
				//Fringe node coupled
				/*				for (int i = 0; i < n_patch_nodes; i++)

				{
					ModelPart::NodesContainerType::iterator iparticle = mrPatchModelPart.NodesBegin() + i;
					Node<3>::Pointer p_patch_node = *(iparticle.base());

					distance = -p_patch_node->FastGetSolutionStepValue(DISTANCE);

					

					if ((distance > 0.4 * overlap_distance) && (distance < 0.6 * overlap_distance))
					{

						typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

						Element::Pointer pElement;

						bool is_found = false;
						is_found = pBinLocator->FindPointOnMesh(p_patch_node->Coordinates(), N, pElement, result_begin, max_results);

						if (is_found == true)
						{

							Geometry<Node<3>> &geom = pElement->GetGeometry();

							std::cout << "Presssure of " << p_patch_node->Id() << "coupled to " << pElement->Id() << std::endl;

							for (int i = 0; i < geom.size(); i++)
							{

								pMpcProcess->AddMasterSlaveRelationWithNodesAndVariable(geom[i], PRESSURE, *p_patch_node, PRESSURE, N[i]);
							}
							IsCoupled = true;
						}

						else
						{
							std::cout << "Cannot find host element for Pressure coupling" << std::endl;
							std::exit(-1);
						}

					} //end of distance check condition
					if(IsCoupled)
					break;
					
				}*/ // end of loop over nodes

				ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin();
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

			} // end of if (type == 0) conditions

			std::cout << "mpcProcess ends" << std::endl;
		}
	}

	void ApplyMpcConstraintConservative(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator, unsigned int type = 1)
	{

		double rtMinvR = 0;
		DofVectorType slaveDofVector;
		double R = 0;

		ApplyMpcConstraint(rBoundaryModelPart, pBinLocator, type);

		std::cout << "Slave nodes are coupled to the nearest element" << std::endl;

		std::vector<VariableComponentType> dofComponentVector = {VELOCITY_X, VELOCITY_Y, VELOCITY_Z};

		// Calculation of Rt*Minv*R and assignment of nodalnormals to the slave dofs

		for (ModelPart::NodesContainerType::iterator inode = rBoundaryModelPart.NodesBegin(); inode != rBoundaryModelPart.NodesEnd(); ++inode)
		{

			double Minode = inode->FastGetSolutionStepValue(NODAL_MASS);

			for (unsigned int i = 0; i < TDim; i++)
			{

				double rIdof = inode->FastGetSolutionStepValue(NORMAL)[i];
				DofType &slaveDOF = inode->GetDof(dofComponentVector[i]);
				pMpcProcess->AddNodalNormalSlaveRelationWithDofs(inode->GetDof(dofComponentVector[i]), rIdof);
				slaveDofVector.push_back(slaveDOF);
				//std::cout << "slaveid,slavekey" << slaveDofVector[i].Id() << "," << slaveDofVector[i].GetVariable().Key() << std::endl;
				rtMinvR += (rIdof * rIdof) / Minode;
				R += rIdof;
			}
		}

		std::cout << "trMinvR calculated : " << rtMinvR << std::endl;
		std::cout << "R calculated : " << R << std::endl;

		SlavePairType slaveDofMap;
		SlavePairType slaveDofMapOther;
		MasterDofWeightMapType MasterDofWeightMap;
		ProcessInfoType info = mrAllModelPart.GetProcessInfo();
		double nodalMass;
		double modifiedWeight;
		double NodalNormalComponent;
		double NodalNormalComponentOther;
		unsigned int slaveNodeId;
		unsigned int slaveDofKey;
		unsigned int slaveNodeIdOther;
		unsigned int slaveDofKeyOther;
		unsigned int masterNodeId;
		unsigned int PartitionId;
		unsigned int masterDofKey;
		unsigned int slaveDofVectorSize = slaveDofVector.size();

		//debug

		/*{
			slaveNodeId = slaveDofVector[0].Id();
			slaveDofKey = slaveDofVector[0].GetVariable().Key();
			slaveDofMap = std::make_pair(slaveNodeId, slaveDofKey);
			MpcDataPointerVectorType mpcDataVector = info.GetValue(MPC_DATA_CONTAINER);
			for (auto mpcData : (*mpcDataVector))
			{

				if (mpcData->IsActive())
				{

					MasterDofWeightMap = mpcData->mDofConstraints[slaveDofMap];
					unsigned int counterbefore = 0;

					for (auto master : MasterDofWeightMap)
					{
						counterbefore++;
						std::tie(masterNodeId, masterDofKey, PartitionId) = master.first;

						std::cout << "master node " << masterNodeId << " coupled to " << slaveNodeId << " weight " << master.second << std::endl;
					}

					std::cout << " Total number of masters " << counterbefore << std::endl;
					std::cout << "############################################################" << std::endl;
				}
			}
		}*/

		MpcDataPointerVectorType mpcDataVector = info.GetValue(MPC_DATA_CONTAINER);
		for (auto mpcData : (*mpcDataVector))
		{

			if (mpcData->IsActive())
			{

				std::cout << "slaveDofVectorSize " << slaveDofVectorSize << std::endl;

				for (unsigned int i = 0; i < slaveDofVectorSize; i++)
				{

					//std::cout<<"Inside i loop nav "<<i <<std::endl;

					slaveNodeId = slaveDofVector[i].Id();
					slaveDofKey = slaveDofVector[i].GetVariable().Key();
					slaveDofMap = std::make_pair(slaveNodeId, slaveDofKey);
					Node<3> &slaveNode = rBoundaryModelPart.Nodes()[slaveNodeId];
					Node<3>::DofsContainerType::iterator it_SlaveDof = slaveNode.GetDofs().find(slaveDofKey);
					nodalMass = slaveNode.FastGetSolutionStepValue(NODAL_MASS);

					for (unsigned int j = 0; j < slaveDofVectorSize; j++)
					{

						//std::cout << "Inside j loop nav " << j << std::endl;

						slaveNodeIdOther = slaveDofVector[j].Id();
						slaveDofKeyOther = slaveDofVector[j].GetVariable().Key();
						slaveDofMapOther = std::make_pair(slaveNodeIdOther, slaveDofKeyOther);

						MasterDofWeightMap = mpcData->mDofConstraints[slaveDofMapOther];
						//std::cout << "slaveId,key " << slaveNodeIdOther << "," << slaveDofKeyOther << std::endl;
						for (auto master : MasterDofWeightMap)

						{

							//std::cout<<"master loop nav "<<std::endl;

							std::tie(masterNodeId, masterDofKey, PartitionId) = master.first;
							Node<3> &masterNode = mrAllModelPart.Nodes()[masterNodeId];
							Node<3>::DofsContainerType::iterator itMasterDof = masterNode.GetDofs().find(masterDofKey);
							NodalNormalComponent = mpcData->mSlaveDofToNodalNormalMap[slaveDofMap];
							NodalNormalComponentOther = mpcData->mSlaveDofToNodalNormalMap[slaveDofMapOther];
							modifiedWeight = -master.second * NodalNormalComponent * NodalNormalComponentOther / (rtMinvR * nodalMass);
							pMpcProcess->AddMasterSlaveRelationWithDofs(*it_SlaveDof, *itMasterDof, modifiedWeight);
							//std::cout << "master ID " << masterNodeId << "," << masterDofKey << " modifiedWeight " << modifiedWeight << std::endl;
						}
					}

					//debug

					MasterDofWeightMap = mpcData->mDofConstraints[slaveDofMap];
					unsigned int counter = 0;
					std::ofstream myfile;
					//myfile.open("example.py");
					//myfile << "l = ";

					for (auto master : MasterDofWeightMap)
					{
						counter++;
						std::tie(masterNodeId, masterDofKey, PartitionId) = master.first;

						//std::cout << "master node " << masterNodeId << "," << masterDofKey << " coupled to " << slaveNodeId << "," << slaveDofKey << " weight " << master.second << std::endl;
						//myfile << masterNodeId << ",";
					}
					//myfile << "\n";
					//myfile << "set([x for x in l if l.count(x) > 1]) \n";
					//myfile.close();

					/*std::cout << " Total number of masters " << counter << std::endl;
					std::cout << "trMinvR calculated : " << rtMinvR << std::endl;
					std::cout << "R calculated : " << R << std::endl;
					std::exit(-1);*/

					//debug
				}
			}
		}
	}

	//Apply Chimera with or without overlap
	void ApplyChimeraUsingMpc(ModelPart &mrPatchBoundaryModelPart, std::string type = "NearestElement")

	{

		for (ModelPart::ElementsContainerType::iterator it = mrAllModelPart.ElementsBegin(); it != mrAllModelPart.ElementsEnd(); ++it)
		{

			it->Set(ACTIVE, true);
		}

		const double epsilon = 1e-12;
		if (overlap_distance < epsilon)
		{
			std::cout << "Overlap distance should be a positive number" << std::endl;
			std::exit(-1);
		}

		if (overlap_distance > epsilon)

		{

			ModelPart::Pointer pHoleModelPart = ModelPart::Pointer(new ModelPart("HoleModelpart"));
			ModelPart::Pointer pHoleBoundaryModelPart = ModelPart::Pointer(new ModelPart("HoleBoundaryModelPart"));
			//ModelPart::Pointer pNewSkinModelPart = ModelPart::Pointer(new ModelPart("NewSkinModelPart"));
			pMpcProcess->SetWeak(true);
			std::cout << "checkpoint nav 1" << std::endl;

			//this->pCalculateDistanceProcess->CalculateSignedDistanceOnModelPart(mrPatchModelPart, mrPatchBoundaryModelPart);
			this->pCalculateDistanceProcess->CalculateSignedDistance(mrBackgroundModelPart, mrPatchBoundaryModelPart);
			//PrintGIDMesh(*pNewSkinModelPart);
			this->pHoleCuttingProcess->CreateHoleAfterDistance(mrBackgroundModelPart, *pHoleModelPart, *pHoleBoundaryModelPart, overlap_distance);

			CalculateNodalAreaAndNodalMass(mrPatchBoundaryModelPart, 1);
			std::cout << "Nodal mass and normal calculated for the patch boundary" << std::endl;
			CalculateNodalAreaAndNodalMass(*pHoleBoundaryModelPart, -1);
			std::cout << "Nodal mass and normal calculated for the hole boundary" << std::endl;

			if (type == "NearestElement")
			{
				ApplyMpcConstraint(mrPatchBoundaryModelPart, pBinLocatorForBackground, 1); //0 for pressure coupling
				std::cout << "Patch boundary coupled with background" << std::endl;
				ApplyMpcConstraint(*pHoleBoundaryModelPart, pBinLocatorForPatch, 1);
				std::cout << "HoleBoundary  coupled with patch" << std::endl;
			}

			else if (type == "Conservative")
			{
				ApplyMpcConstraintConservative(mrPatchBoundaryModelPart, pBinLocatorForBackground, 0); //0 for pressure coupling
				std::cout << "Patch boundary coupled with background using conservative approach" << std::endl;
				ApplyMpcConstraintConservative(*pHoleBoundaryModelPart, pBinLocatorForPatch, 1);
				std::cout << "HoleBoundary  coupled with patch using conservative approach" << std::endl;
			}
		}

		else
		{

			std::cout << "Applying Strong MPC " << std::endl;

			ModelPart::Pointer pHoleModelPart = ModelPart::Pointer(new ModelPart("HoleModelpart"));
			ModelPart::Pointer pHoleBoundaryModelPart = ModelPart::Pointer(new ModelPart("HoleBoundaryModelPart"));
			//ModelPart::Pointer pNewSkinModelPart = ModelPart::Pointer(new ModelPart("NewSkinModelPart"));

			//this->pCalculateDistanceProcess->ExtractDistance(mrPatchModelPart, mrBackgroundModelPart, mrPatchBoundaryModelPart);
			this->pCalculateDistanceProcess->CalculateSignedDistance(mrBackgroundModelPart, mrPatchBoundaryModelPart);

			this->pHoleCuttingProcess->CreateHoleAfterDistance(mrBackgroundModelPart, *pHoleModelPart, *pHoleBoundaryModelPart, epsilon);

			if (type == "NearestElement")
				ApplyMpcConstraint(mrPatchBoundaryModelPart, pBinLocatorForBackground, 1);

			else if (type == "Conservative")
				ApplyMpcConstraintConservative(mrPatchBoundaryModelPart, pBinLocatorForBackground, 1);
		}
	}

	void SetOverlapDistance(double distance)
	{

		this->overlap_distance = distance;
	}

	void CalculateNodalAreaAndNodalMass(ModelPart &rBoundaryModelPart, int sign)
	{
		KRATOS_TRY
		ConditionsArrayType &rConditions = rBoundaryModelPart.Conditions();
		//resetting the normals and calculating centre point

		array_1d<double, 3> zero;
		array_1d<double, 3> centre;
		unsigned int n_nodes = rBoundaryModelPart.Nodes().size();

		zero[0] = 0.0;
		zero[1] = 0.0;
		zero[2] = 0.0;

		centre[0] = 0.0;
		centre[1] = 0.0;
		centre[2] = 0.0;

		for (ConditionsArrayType::iterator it = rConditions.begin();
			 it != rConditions.end(); it++)
		{
			Element::GeometryType &rNodes = it->GetGeometry();
			for (unsigned int in = 0; in < rNodes.size(); in++)
			{
				noalias((rNodes[in]).GetSolutionStepValue(NORMAL)) = zero;
			}
		}

		for (ModelPart::NodesContainerType::iterator inode = rBoundaryModelPart.NodesBegin(); inode != rBoundaryModelPart.NodesEnd(); ++inode)
		{

			centre += inode->Coordinates();
		}

		centre = centre / n_nodes;
		std::cout << "Centre " << centre[0] << " " << centre[1] << " " << centre[2] << std::endl;

		//calculating the normals and storing on the conditions
		array_1d<double, 3> An;
		if (TDim == 2)
		{
			for (ConditionsArrayType::iterator it = rConditions.begin();
				 it != rConditions.end(); it++)
			{
				if (it->GetGeometry().PointsNumber() == 2)
					CalculateNormal2D(it, An, centre, sign);
			}
		}
		else if (TDim == 3)
		{
			array_1d<double, 3> v1;
			array_1d<double, 3> v2;
			for (ConditionsArrayType::iterator it = rConditions.begin();
				 it != rConditions.end(); it++)
			{
				//calculate the normal on the given condition
				if (it->GetGeometry().PointsNumber() == 3)
					CalculateNormal3D(it, An, v1, v2, centre, sign);
			}
		}

		//adding the normals to the nodes
		for (ConditionsArrayType::iterator it = rConditions.begin();
			 it != rConditions.end(); it++)
		{
			Geometry<Node<3>> &pGeometry = (it)->GetGeometry();
			double coeff = 1.00 / pGeometry.size();
			const array_1d<double, 3> &normal = it->GetValue(NORMAL);
			double nodal_mass = MathUtils<double>::Norm3(normal);

			for (unsigned int i = 0; i < pGeometry.size(); i++)
			{
				noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * normal;
				pGeometry[i].FastGetSolutionStepValue(NODAL_MASS) += coeff * nodal_mass;
			}
		}

		KRATOS_CATCH("")
	}

	void CalculateNormal2D(ConditionsArrayType::iterator it, array_1d<double, 3> &An, array_1d<double, 3> &centre, int sign)
	{
		Geometry<Node<3>> &pGeometry = (it)->GetGeometry();
		array_1d<double, 3> rVector;

		An[0] = pGeometry[1].Y() - pGeometry[0].Y();
		An[1] = -(pGeometry[1].X() - pGeometry[0].X());
		An[2] = 0.00;

		rVector[0] = centre[0] - pGeometry[0].X();
		rVector[1] = centre[1] - pGeometry[0].Y();
		rVector[2] = 0.00;

		array_1d<double, 3> &normal = (it)->GetValue(NORMAL);
		noalias(normal) = An;

		if ((MathUtils<double>::Dot(An, rVector) > 0))
			normal = -1 * normal;

		normal = normal * sign;

		// 				(it)->SetValue(NORMAL,An);
	}

	void CalculateNormal3D(ConditionsArrayType::iterator it, array_1d<double, 3> &An,
						   array_1d<double, 3> &v1, array_1d<double, 3> &v2, array_1d<double, 3> &centre, int sign)
	{
		Geometry<Node<3>> &pGeometry = (it)->GetGeometry();
		array_1d<double, 3> rVector;

		v1[0] = pGeometry[1].X() - pGeometry[0].X();
		v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
		v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

		v2[0] = pGeometry[2].X() - pGeometry[0].X();
		v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
		v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

		rVector[0] = centre[0] - pGeometry[0].X();
		rVector[1] = centre[1] - pGeometry[0].Y();
		rVector[2] = centre[2] - pGeometry[0].Z();

		MathUtils<double>::CrossProduct(An, v1, v2);
		An *= 0.5;

		array_1d<double, 3> &normal = (it)->GetValue(NORMAL);
		noalias(normal) = An;

		if ((MathUtils<double>::Dot(An, rVector) > 0))
			normal = -1 * normal;

		normal = normal * sign;
		// 				noalias((it)->GetValue(NORMAL)) = An;
	}

	void PrintGIDMesh(ModelPart& rmodel_part)
	{
		std::ofstream myfile;
        myfile.open (rmodel_part.Name()+".post.msh");
        myfile << "MESH \"leaves\" dimension 2 ElemType Line Nnode 2" << std::endl;
        myfile << "# color 96 96 96" << std::endl;
        myfile << "Coordinates" << std::endl;
        myfile << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;
        
		

		for (unsigned int i = 0; i < rmodel_part.Nodes().size(); i++)
		{
			ModelPart::NodesContainerType::iterator iparticle = rmodel_part.NodesBegin() + i;
			Node<3>::Pointer p_node = *(iparticle.base());
			myfile << p_node->Id()<< "  " << p_node->Coordinates()[0] << "  " << p_node->Coordinates()[1]<< "  " << p_node->Coordinates()[2] << std::endl;

		}

		myfile << "end coordinates" << std::endl;
        myfile << "elements" << std::endl;
        myfile << "# element node_1 node_2 material_number" << std::endl;

		for (ConditionsArrayType::iterator it = rmodel_part.Conditions().begin();
			 it != rmodel_part.Conditions().end(); it++)
		{
			
			myfile << it->Id() << "  ";
			for (unsigned int i =0; i < it->GetGeometry().PointsNumber(); i++)
			myfile << (it->GetGeometry()[i]).Id() << "  ";

			myfile << std::endl;

		}

		myfile << "end elements" << std::endl;




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
