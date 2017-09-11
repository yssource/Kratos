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

#if !defined(TEST_MAPPER_H_INCLUDED)
#define TEST_MAPPER_H_INCLUDED

// System includes

#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

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
#include "includes/kratos_parameters.h"

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

class TestMapperProcess
{
  public:
	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomExtractVariablesProcess
	KRATOS_CLASS_POINTER_DEFINITION(TestMapperProcess);

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
	TestMapperProcess(ModelPart &AllModelPart, ModelPart &BackgroundModelPart, ModelPart &PatchModelPart, double distance = 1e-12) : mrAllModelPart(AllModelPart), mrBackgroundModelPart(BackgroundModelPart), mrPatchModelPart(PatchModelPart), overlap_distance(distance)

	{
		this->pBinLocatorForBackground = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(mrBackgroundModelPart));

		//this->pMpcProcess = ApplyMultipointConstraintsProcess::Pointer(new ApplyMultipointConstraintsProcess(mrAllModelPart));
		this->pMpcProcess = NULL; //mapping at the hole boundary
		this->pCalculateDistanceProcess = typename CustomCalculateSignedDistanceProcess<TDim>::Pointer(new CustomCalculateSignedDistanceProcess<TDim>());
	}

	/// Destructor.
	virtual ~TestMapperProcess()
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
			;
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

			// For one ndoe pressure coupling
			/*{
				ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin();
				Node<3>::Pointer p_boundary_node = *(iparticle.base());

				typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

				Element::Pointer pElement;

				bool is_found = false;
				is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

				if (is_found == true)
				{

					Geometry<Node<3>> &geom = pElement->GetGeometry();

					{

						for (int i = 0; i < geom.size(); i++)
						{

							pMpcProcess->AddMasterSlaveRelationWithNodesAndVariable(geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);
						}
					}
				}
			}
*/
			std::cout << "mpcProcess ends" << std::endl;
		}
	}

	void ApplyMpcConstraintConservative(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator, unsigned int type = 1)
	{

		double rtMinvR = 0;

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
				pMpcProcess->AddNodalNormalSlaveRelationWithDofs(inode->GetDof(dofComponentVector[i]), rIdof);

				//std::cout << "slaveid,slavekey" << slaveDofVector[i].Id() << "," << slaveDofVector[i].GetVariable().Key() << std::endl;
				rtMinvR += (rIdof * rIdof) / Minode;
			}

			pMpcProcess->AddNodalNormalSlaveRelationWithDofs(inode->GetDof(PRESSURE), 0);

			pMpcProcess->SetRtMinvR(rtMinvR);
		}

		/*// One node pressure coupling
		ModelPart::NodesContainerType::iterator inode = rBoundaryModelPart.NodesBegin();
		pMpcProcess->AddNodalNormalSlaveRelationWithDofs(inode->GetDof(PRESSURE), 0);*/
	}

	/*void Initialize(ModelPart& rPatchBoundaryModelPart, std::string type = "NearestElement")
	{

		std::ofstream flux_file;
		flux_file.open("background_flux.csv");
		flux_file << "Time,massflux" << std::endl;

		this->pCalculateDistanceProcess->CalculateSignedDistance(mrBackgroundModelPart, rPatchBoundaryModelPart);
		
		CalculateNodalAreaAndNodalMass(rPatchBoundaryModelPart, 1);
		std::cout << "Nodal mass and normal calculated for the patch boundary" << std::endl;

		if (type == "NearestElement")
		{
			ApplyMpcConstraint(rPatchBoundaryModelPart, pBinLocatorForBackground, 1); //0 for pressure coupling
			std::cout << "Patch boundary coupled with background" << std::endl;
		}

		else if (type == "Conservative")
		{
			ApplyMpcConstraintConservative(rPatchBoundaryModelPart, pBinLocatorForBackground, 1); //0 for pressure coupling
			std::cout << "Patch boundary coupled with background using conservative approach" << std::endl;
		}

	
	}*/

	//Apply Chimera with or without overlap
	/*void Interpolate2d(ModelPart &rPatchBoundaryModelPart, std::string type = "NearestElement", double t = 0)

	{

	

		//Calculate R

		ProcessInfo &CurrentProcessInfo = mrAllModelPart.GetProcessInfo();
		MpcDataPointerVectorType& mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

		for (auto mpcData : (*mpcDataVector))
		{

			if (mpcData->IsActive())

			{
				double r = 0;

				//#####################for Creating skin model part#####################################################
				unsigned int id_node = 1;
				unsigned int id_condition = 1;
				ModelPart::Pointer pSkinModelPart = ModelPart::Pointer(new ModelPart("SkinModelPart"));
				//#####################for Creating skin model part#####################################################
				for (ModelPart::ElementsContainerType::iterator i_fluid_element = mrBackgroundModelPart.ElementsBegin(); i_fluid_element != mrBackgroundModelPart.ElementsEnd(); i_fluid_element++)
				{
					double NumberOfPostiveDistance = 0;
					bool &is_split = i_fluid_element->GetValue(SPLIT_ELEMENT);
					Geometry<Node<3>> &geom = i_fluid_element->GetGeometry();
					for (unsigned int i = 0; i < geom.size(); i++)
					{

						double distance = geom[i].FastGetSolutionStepValue(DISTANCE);

						if (distance > 0)
							NumberOfPostiveDistance++;
					}

					if (NumberOfPostiveDistance == geom.size())
					{
						is_split = false;
					}
				}

				for (ModelPart::ElementsContainerType::iterator i_fluid_element = mrBackgroundModelPart.ElementsBegin(); i_fluid_element != mrBackgroundModelPart.ElementsEnd(); i_fluid_element++)
				{
					bool is_split = i_fluid_element->GetValue(SPLIT_ELEMENT);

					if (is_split == true)
					{

						Geometry<Node<3>> &geom = i_fluid_element->GetGeometry();
						unsigned int numberOfPoints = geom.size();
						std::vector<double> distances(numberOfPoints, 0.0);

						for (unsigned int i = 0; i < numberOfPoints; i++)
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);

						// generate the points on the edges at the zero of the distance function
						std::vector<Point<3>> edge_points;
						std::vector<double> vectorOfVelx;
						std::vector<double> vectorOfVely;
						edge_points.reserve(2);
						array_1d<double, 3> rVector;
						double d_nodei, d_nodej;
						array_1d<double, 3> An;

						// loop over all 3 edges of the triangle
						for (unsigned int i = 0; i < 2; i++)
						{
							for (unsigned int j = i + 1; j < 3; j++) // go through the edges 01, 02, 12
							{
								double di = distances[i];
								double dj = distances[j];

								if (di * dj < 0) //edge is cut
								{
									// generate point on edge by linear interpolation
									double Ni = fabs(dj) / (fabs(di) + fabs(dj));
									double Nj = 1.0 - Ni;
									Point<3> edge_point(Ni * geom[i] + Nj * geom[j]);
									edge_points.push_back(edge_point);
									double velx = geom[i].FastGetSolutionStepValue(VELOCITY_X) * Ni + geom[j].FastGetSolutionStepValue(VELOCITY_X) * Nj;
									double vely = geom[i].FastGetSolutionStepValue(VELOCITY_Y) * Ni + geom[j].FastGetSolutionStepValue(VELOCITY_Y) * Nj;

									vectorOfVelx.push_back(velx);
									vectorOfVely.push_back(vely);
									d_nodej = geom[j].FastGetSolutionStepValue(DISTANCE);
									d_nodei = geom[i].FastGetSolutionStepValue(DISTANCE);

									//rVector to get the outward orientation of normal
									rVector[0] = (d_nodej - d_nodei) * (geom[j].X() - geom[i].X());
									rVector[1] = (d_nodej - d_nodei) * (geom[j].Y() - geom[i].Y());
									rVector[2] = 0.0;
								}
							}
						}
						// Normal calculation
						An[0] = edge_points[1].Y() - edge_points[0].Y();
						An[1] = -(edge_points[1].X() - edge_points[0].X());
						An[2] = 0.00;

						if ((MathUtils<double>::Dot(An, rVector) < 0))
							An = -1 * An;

						array_1d<double, 3> &normal = i_fluid_element->GetValue(NORMAL);
						normal = An;
						//###################################FOR VISUALSING SKIN MODEL PART##################################################################################
						if (edge_points.size() == 2)
						{
							// ######## ADDING NEW NODE #########

							Node<3>::Pointer pnode1 = pSkinModelPart->CreateNewNode(id_node++, edge_points[0].X(), edge_points[0].Y(), edge_points[0].Z());
							Node<3>::Pointer pnode2 = pSkinModelPart->CreateNewNode(id_node++, edge_points[1].X(), edge_points[1].Y(), edge_points[1].Z());

							// ######## ADDING NEW CONDITION #########
							//form a triangle
							Line2D2<Node<3>> line(pnode1, pnode2);

							Condition const &rReferenceCondition = KratosComponents<Condition>::Get("Condition2D");
							Properties::Pointer properties = pSkinModelPart->rProperties()(0);
							Condition::Pointer p_condition = rReferenceCondition.Create(id_condition++, line, properties);

							pSkinModelPart->Conditions().push_back(p_condition);
						}
						PrintGIDMesh(*pSkinModelPart);
						//###################################FOR VISUALSING SKIN MODEL PART##################################################################################

						if (edge_points.size() == 2)
						{

							for (unsigned int i = 0; i < 2; i++)
							{

								r += vectorOfVelx[i] * An[0] * 0.5;
								r += vectorOfVely[i] * An[1] * 0.5;
							}
						}

					} // if split = true

				} // End of element loop

				// printing flux from fluid model part

				std::ofstream flux_file;
				
				flux_file.open("background_flux.csv", std::ios_base::app);

				flux_file << t << "," << r << std::endl;

				flux_file.close();
				// Nearest Element Mapping
				for (auto slaveMasterDofMap : mpcData->mDofConstraints)
				{
					SlavePairType slaveDofMap = slaveMasterDofMap.first;
					MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
					unsigned int slaveNodeId = slaveDofMap.first;
					unsigned int slaveDofKey = slaveDofMap.second;
					Node<3> &slaveNode = mrAllModelPart.Nodes()[slaveNodeId];
					Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
					double &slaveDofValue = idof->GetSolutionStepValue();
					slaveDofValue = 0.0;
					
					for (auto masterDofMapElem : masterDofMap)
					{
						unsigned int masterNodeId;
						unsigned int masterDofKey;
						unsigned int PartitionId;
						double weight = masterDofMapElem.second;
						std::tie(masterNodeId, masterDofKey, PartitionId) = masterDofMapElem.first;
						Node<3> &masterNode = mrBackgroundModelPart.Nodes()[masterNodeId];
						Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);
						slaveDofValue += itMaster->GetSolutionStepValue() * weight;
						
					} // masterDofMapElem loop
					

				} // slaveMasterDofMap loop

				std::cout << "Finished mapping by Nearest Element" << std::endl;
				// For Conservative

				if (type == "Conservative")
				{
					double nodalMass;
					unsigned int slaveNodeId;
					unsigned int slaveNodeIdOther;
					unsigned int slaveDofKey;
					unsigned int slaveDofKeyOther;
					double slaveDofValueOther;
					SlavePairType slaveDofMap;
					SlavePairType slaveDofMapOther;
					double RtMinvR = mpcData->RtMinvR;
					double NodalNormalComponent;
					double NodalNormalComponentOther;
					std::cout << " RtMinvR " << RtMinvR << std::endl;
					std::vector<double> VectorOfconstants;
					unsigned int slaveIndex = 0;

					for (auto slaveMasterDofMap : mpcData->mDofConstraints)
					{
						slaveDofMap = slaveMasterDofMap.first;
						slaveNodeId = slaveDofMap.first;
						slaveDofKey = slaveDofMap.second;
						Node<3> &slaveNode = mrAllModelPart.Nodes()[slaveNodeId];
						Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
						nodalMass = slaveNode.FastGetSolutionStepValue(NODAL_MASS);
						NodalNormalComponent = mpcData->mSlaveDofToNodalNormalMap[slaveDofMap];
						VectorOfconstants.push_back(0.0);
						for (auto slaveMasterDofMapOther : mpcData->mDofConstraints)
						{

							slaveDofMapOther = slaveMasterDofMapOther.first;
							slaveNodeIdOther = slaveDofMapOther.first;
							slaveDofKeyOther = slaveDofMapOther.second;
							Node<3> &slaveNodeOther = mrAllModelPart.Nodes()[slaveNodeIdOther];
							Node<3>::DofsContainerType::iterator idofOther = slaveNodeOther.GetDofs().find(slaveDofKeyOther);
							slaveDofValueOther = idofOther->GetSolutionStepValue();
							NodalNormalComponentOther = mpcData->mSlaveDofToNodalNormalMap[slaveDofMapOther];
							VectorOfconstants[slaveIndex] -= ((NodalNormalComponent * NodalNormalComponentOther) / (nodalMass * RtMinvR)) * slaveDofValueOther; // correction for zero flux

							//debug
							//if((slaveDofKey == 972)&&(slaveDofKeyOther == 972))
								//std::cout << "Correction for 0 flux on slave: " << slaveNodeId << " from " << slaveNodeIdOther <<" :: "<<((NodalNormalComponent * NodalNormalComponentOther) / (nodalMass * RtMinvR)) <<std::endl;

						} // slaveMasterDofMap loop

						VectorOfconstants[slaveIndex] += ((NodalNormalComponent * r) / (nodalMass * RtMinvR));
						slaveIndex++;
						//debug
						//if((slaveDofKey == 972))
						//std::cout << "Correction for finite flux :: " << ((NodalNormalComponent * r) / (nodalMass * RtMinvR)) << std::endl;

					} // slaveMasterDofMapOther loop

					slaveIndex = 0;

					//Applying correction in the slaveDofValue
					for (auto slaveMasterDofMap : mpcData->mDofConstraints)
					{

						slaveDofMap = slaveMasterDofMap.first;
						slaveNodeId = slaveDofMap.first;
						slaveDofKey = slaveDofMap.second;
						Node<3> &slaveNode = mrAllModelPart.Nodes()[slaveNodeId];
						Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
						double &slaveDofValue = idof->GetSolutionStepValue();

						slaveDofValue += VectorOfconstants[slaveIndex];
						slaveIndex++;

					} // slaveMasterDofMap loop

					std::cout << "Conservative correction to the velocity field applied" << std::endl;

				} // if type == "Conservative"

			} // mpcData->IsActive()

		} // mpcData vector

		//mpcDataVector = NULL;

	} // End of function*/

	void Initialize()
	{

		std::ofstream flux_file;
		flux_file.open("background_flux.csv");
		flux_file << "Time,massflux" << std::endl;
	}

	void Interpolate2d(ModelPart &rPatchBoundaryModelPart, std::string type = "NearestElement", double t = 0)

	{
		ModelPart::Pointer pHoleModelPart = ModelPart::Pointer(new ModelPart("HoleModelpart"));
		ModelPart::Pointer pHoleBoundaryModelPart = ModelPart::Pointer(new ModelPart("HoleBoundaryModelpart"));
		this->pMpcProcess = ApplyMultipointConstraintsProcess::Pointer(new ApplyMultipointConstraintsProcess(mrAllModelPart));

		for (ModelPart::ElementsContainerType::iterator i_elem = mrAllModelPart.ElementsBegin(); i_elem != mrAllModelPart.ElementsBegin(); i_elem++)
		{

			i_elem->Set(ACTIVE, true);
		}

		this->pCalculateDistanceProcess->CalculateSignedDistance(mrPatchModelPart, rPatchBoundaryModelPart);
		this->pHoleCuttingProcess->CreateHoleAfterDistance(mrPatchModelPart, *pHoleModelPart, *pHoleBoundaryModelPart, overlap_distance);
		this->pCalculateDistanceProcess->CalculateSignedDistance(mrBackgroundModelPart, *pHoleBoundaryModelPart);

		CalculateNodalAreaAndNodalMass(*pHoleBoundaryModelPart, 1);

		std::cout << "Nodal mass and normal calculated for the patch boundary" << std::endl;
	
		if (type == "NearestElement")
		{
			ApplyMpcConstraint(*pHoleBoundaryModelPart, pBinLocatorForBackground, 1); //0 for pressure coupling
			std::cout << "Patch boundary coupled with background" << std::endl;
		}

		else if (type == "Conservative")
		{
			ApplyMpcConstraintConservative(*pHoleBoundaryModelPart, pBinLocatorForBackground, 1); //0 for pressure coupling
			std::cout << "Patch boundary coupled with background using conservative approach" << std::endl;
		}

		for (ModelPart::NodesContainerType::iterator i_node = pHoleBoundaryModelPart->NodesBegin(); i_node != pHoleBoundaryModelPart->NodesEnd(); i_node++)
		{

			i_node->Fix(VELOCITY_X);
			i_node->Fix(VELOCITY_Y);
			i_node->Fix(PRESSURE);
		}

		//Calculate R

		ProcessInfo &CurrentProcessInfo = mrAllModelPart.GetProcessInfo();
		MpcDataPointerVectorType &mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

		for (auto mpcData : (*mpcDataVector))
		{

			if (mpcData->IsActive())

			{
				double r = 0;

				//#####################for Creating skin model part#####################################################
				unsigned int id_node = 1;
				unsigned int id_condition = 1;
				ModelPart::Pointer pSkinModelPart = ModelPart::Pointer(new ModelPart("SkinModelPart"));
				//#####################for Creating skin model part#####################################################
				for (ModelPart::ElementsContainerType::iterator i_fluid_element = mrBackgroundModelPart.ElementsBegin(); i_fluid_element != mrBackgroundModelPart.ElementsEnd(); i_fluid_element++)
				{
					double NumberOfPostiveDistance = 0;
					bool &is_split = i_fluid_element->GetValue(SPLIT_ELEMENT);
					Geometry<Node<3>> &geom = i_fluid_element->GetGeometry();
					for (unsigned int i = 0; i < geom.size(); i++)
					{

						double distance = geom[i].FastGetSolutionStepValue(DISTANCE);

						if (distance > 0)
							NumberOfPostiveDistance++;
					}

					if (NumberOfPostiveDistance == geom.size())
					{
						is_split = false;
					}
				}

				for (ModelPart::ElementsContainerType::iterator i_fluid_element = mrBackgroundModelPart.ElementsBegin(); i_fluid_element != mrBackgroundModelPart.ElementsEnd(); i_fluid_element++)
				{
					bool is_split = i_fluid_element->GetValue(SPLIT_ELEMENT);

					if (is_split == true)
					{

						Geometry<Node<3>> &geom = i_fluid_element->GetGeometry();
						unsigned int numberOfPoints = geom.size();
						std::vector<double> distances(numberOfPoints, 0.0);

						for (unsigned int i = 0; i < numberOfPoints; i++)
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);

						// generate the points on the edges at the zero of the distance function
						std::vector<Point<3>> edge_points;
						std::vector<double> vectorOfVelx;
						std::vector<double> vectorOfVely;
						edge_points.reserve(2);
						array_1d<double, 3> rVector;
						double d_nodei, d_nodej;
						array_1d<double, 3> An;

						// loop over all 3 edges of the triangle
						for (unsigned int i = 0; i < 2; i++)
						{
							for (unsigned int j = i + 1; j < 3; j++) // go through the edges 01, 02, 12
							{
								double di = distances[i];
								double dj = distances[j];

								if (di * dj < 0) //edge is cut
								{
									// generate point on edge by linear interpolation
									double Ni = fabs(dj) / (fabs(di) + fabs(dj));
									double Nj = 1.0 - Ni;
									Point<3> edge_point(Ni * geom[i] + Nj * geom[j]);
									edge_points.push_back(edge_point);
									double velx = geom[i].FastGetSolutionStepValue(VELOCITY_X) * Ni + geom[j].FastGetSolutionStepValue(VELOCITY_X) * Nj;
									double vely = geom[i].FastGetSolutionStepValue(VELOCITY_Y) * Ni + geom[j].FastGetSolutionStepValue(VELOCITY_Y) * Nj;

									vectorOfVelx.push_back(velx);
									vectorOfVely.push_back(vely);
									d_nodej = geom[j].FastGetSolutionStepValue(DISTANCE);
									d_nodei = geom[i].FastGetSolutionStepValue(DISTANCE);

									//rVector to get the outward orientation of normal
									rVector[0] = (d_nodej - d_nodei) * (geom[j].X() - geom[i].X());
									rVector[1] = (d_nodej - d_nodei) * (geom[j].Y() - geom[i].Y());
									rVector[2] = 0.0;
								}
							}
						}
						// Normal calculation
						An[0] = edge_points[1].Y() - edge_points[0].Y();
						An[1] = -(edge_points[1].X() - edge_points[0].X());
						An[2] = 0.00;

						if ((MathUtils<double>::Dot(An, rVector) < 0))
							An = -1 * An;

						array_1d<double, 3> &normal = i_fluid_element->GetValue(NORMAL);
						normal = An;
						//###################################FOR VISUALSING SKIN MODEL PART##################################################################################
						if (edge_points.size() == 2)
						{
							// ######## ADDING NEW NODE #########

							Node<3>::Pointer pnode1 = pSkinModelPart->CreateNewNode(id_node++, edge_points[0].X(), edge_points[0].Y(), edge_points[0].Z());
							Node<3>::Pointer pnode2 = pSkinModelPart->CreateNewNode(id_node++, edge_points[1].X(), edge_points[1].Y(), edge_points[1].Z());

							// ######## ADDING NEW CONDITION #########
							//form a triangle
							Line2D2<Node<3>> line(pnode1, pnode2);

							Condition const &rReferenceCondition = KratosComponents<Condition>::Get("Condition2D");
							Properties::Pointer properties = pSkinModelPart->rProperties()(0);
							Condition::Pointer p_condition = rReferenceCondition.Create(id_condition++, line, properties);

							pSkinModelPart->Conditions().push_back(p_condition);
						}
						PrintGIDMesh(*pSkinModelPart);
						//###################################FOR VISUALSING SKIN MODEL PART##################################################################################

						if (edge_points.size() == 2)
						{

							for (unsigned int i = 0; i < 2; i++)
							{

								r += vectorOfVelx[i] * An[0] * 0.5;
								r += vectorOfVely[i] * An[1] * 0.5;
							}
						}

					} // if split = true

				} // End of element loop

				// printing flux from fluid model part

				std::ofstream flux_file;

				flux_file.open("background_flux.csv", std::ios_base::app);

				flux_file << t << "," << r << std::endl;

				flux_file.close();
				// Nearest Element Mapping
				for (auto slaveMasterDofMap : mpcData->mDofConstraints)
				{
					SlavePairType slaveDofMap = slaveMasterDofMap.first;
					MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
					unsigned int slaveNodeId = slaveDofMap.first;
					unsigned int slaveDofKey = slaveDofMap.second;
					Node<3> &slaveNode = mrAllModelPart.Nodes()[slaveNodeId];
					Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
					double &slaveDofValue = idof->GetSolutionStepValue();
					slaveDofValue = 0.0;

					for (auto masterDofMapElem : masterDofMap)
					{
						unsigned int masterNodeId;
						unsigned int masterDofKey;
						unsigned int PartitionId;
						double weight = masterDofMapElem.second;
						std::tie(masterNodeId, masterDofKey, PartitionId) = masterDofMapElem.first;
						Node<3> &masterNode = mrBackgroundModelPart.Nodes()[masterNodeId];
						Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);
						slaveDofValue += itMaster->GetSolutionStepValue() * weight;

					} // masterDofMapElem loop

				} // slaveMasterDofMap loop

				std::cout << "Finished mapping by Nearest Element" << std::endl;
				// For Conservative

				if (type == "Conservative")
				{
					double nodalMass;
					unsigned int slaveNodeId;
					unsigned int slaveNodeIdOther;
					unsigned int slaveDofKey;
					unsigned int slaveDofKeyOther;
					double slaveDofValueOther;
					SlavePairType slaveDofMap;
					SlavePairType slaveDofMapOther;
					double RtMinvR = mpcData->RtMinvR;
					double NodalNormalComponent;
					double NodalNormalComponentOther;
					std::cout << " RtMinvR " << RtMinvR << std::endl;
					std::vector<double> VectorOfconstants;
					unsigned int slaveIndex = 0;

					for (auto slaveMasterDofMap : mpcData->mDofConstraints)
					{
						slaveDofMap = slaveMasterDofMap.first;
						slaveNodeId = slaveDofMap.first;
						slaveDofKey = slaveDofMap.second;
						Node<3> &slaveNode = mrAllModelPart.Nodes()[slaveNodeId];
						Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
						nodalMass = slaveNode.FastGetSolutionStepValue(NODAL_MASS);
						NodalNormalComponent = mpcData->mSlaveDofToNodalNormalMap[slaveDofMap];
						VectorOfconstants.push_back(0.0);
						for (auto slaveMasterDofMapOther : mpcData->mDofConstraints)
						{

							slaveDofMapOther = slaveMasterDofMapOther.first;
							slaveNodeIdOther = slaveDofMapOther.first;
							slaveDofKeyOther = slaveDofMapOther.second;
							Node<3> &slaveNodeOther = mrAllModelPart.Nodes()[slaveNodeIdOther];
							Node<3>::DofsContainerType::iterator idofOther = slaveNodeOther.GetDofs().find(slaveDofKeyOther);
							slaveDofValueOther = idofOther->GetSolutionStepValue();
							NodalNormalComponentOther = mpcData->mSlaveDofToNodalNormalMap[slaveDofMapOther];
							VectorOfconstants[slaveIndex] -= ((NodalNormalComponent * NodalNormalComponentOther) / (nodalMass * RtMinvR)) * slaveDofValueOther; // correction for zero flux

							//debug
							//if((slaveDofKey == 972)&&(slaveDofKeyOther == 972))
							//std::cout << "Correction for 0 flux on slave: " << slaveNodeId << " from " << slaveNodeIdOther <<" :: "<<((NodalNormalComponent * NodalNormalComponentOther) / (nodalMass * RtMinvR)) <<std::endl;

						} // slaveMasterDofMap loop

						VectorOfconstants[slaveIndex] += ((NodalNormalComponent * r) / (nodalMass * RtMinvR));
						slaveIndex++;
						//debug
						//if((slaveDofKey == 972))
						//std::cout << "Correction for finite flux :: " << ((NodalNormalComponent * r) / (nodalMass * RtMinvR)) << std::endl;

					} // slaveMasterDofMapOther loop

					slaveIndex = 0;

					//Applying correction in the slaveDofValue
					for (auto slaveMasterDofMap : mpcData->mDofConstraints)
					{

						slaveDofMap = slaveMasterDofMap.first;
						slaveNodeId = slaveDofMap.first;
						slaveDofKey = slaveDofMap.second;
						Node<3> &slaveNode = mrAllModelPart.Nodes()[slaveNodeId];
						Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
						double &slaveDofValue = idof->GetSolutionStepValue();

						slaveDofValue += VectorOfconstants[slaveIndex];
						slaveIndex++;

					} // slaveMasterDofMap loop

					std::cout << "Conservative correction to the velocity field applied" << std::endl;

				} // if type == "Conservative"

			} // mpcData->IsActive()

		} // mpcData vector

		mpcDataVector = NULL;

	} // End of function*/

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
				rNodes[in].GetSolutionStepValue(NODAL_MASS) = 0.0;
			}
		}

		for (ModelPart::NodesContainerType::iterator inode = rBoundaryModelPart.NodesBegin(); inode != rBoundaryModelPart.NodesEnd(); ++inode)
		{

			centre += inode->Coordinates();
		}

		centre = centre / n_nodes;

		/*//hardcoded for Test
		centre[0] = 0.5;
		centre[1] = 0.5;*/
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

	void PrintGIDMesh(ModelPart &rmodel_part)
	{

		std::ofstream myfile;
		myfile.open(rmodel_part.Name() + ".post.msh");
		myfile << "MESH \"leaves\" dimension 2 ElemType Line Nnode 2" << std::endl;
		myfile << "# color 96 96 96" << std::endl;
		myfile << "Coordinates" << std::endl;
		myfile << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

		for (unsigned int i = 0; i < rmodel_part.Nodes().size(); i++)
		{
			ModelPart::NodesContainerType::iterator iparticle = rmodel_part.NodesBegin() + i;
			Node<3>::Pointer p_node = *(iparticle.base());
			myfile << p_node->Id() << "  " << p_node->Coordinates()[0] << "  " << p_node->Coordinates()[1] << "  " << p_node->Coordinates()[2] << std::endl;
		}

		myfile << "end coordinates" << std::endl;
		myfile << "elements" << std::endl;
		myfile << "# element node_1 node_2 material_number" << std::endl;

		for (ConditionsArrayType::iterator it = rmodel_part.Conditions().begin();
			 it != rmodel_part.Conditions().end(); it++)
		{

			myfile << it->Id() << "  ";
			for (unsigned int i = 0; i < it->GetGeometry().PointsNumber(); i++)
				myfile << (it->GetGeometry()[i]).Id() << "  ";

			myfile << std::endl;
		}

		myfile << "end elements" << std::endl;
	}

	virtual std::string Info() const
	{
		return "TestMapperProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "TestMapperProcess";
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
	TestMapperProcess &operator=(TestMapperProcess const &rOther);

	/// Copy constructor.
	//CustomExtractVariablesProcess(CustomExtractVariablesProcess const& rOther);

	///@}
}; // Class CustomExtractVariablesProcess

} // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
