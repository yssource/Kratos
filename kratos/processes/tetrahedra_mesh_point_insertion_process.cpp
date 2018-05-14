//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/exception.h"
#include "includes/kratos_flags.h"
#include "processes/tetrahedra_mesh_point_insertion_process.h"
#include "processes/measure_mesh_quality_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "modeler/tetrahedra_ball.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos
{

	TetrahedraMeshPointInsertionProcess::TetrahedraMeshPointInsertionProcess(ModelPart& rModelPart, double AptQuality, std::size_t IterationsNumber)
		:TetrahedraMeshWorstElementSmoothingProcess(rModelPart, AptQuality, IterationsNumber) {
	}

	TetrahedraMeshPointInsertionProcess::~TetrahedraMeshPointInsertionProcess() {
	}

	/// Turn back information as a string.
	std::string TetrahedraMeshPointInsertionProcess::Info() const
	{
		return "TetrahedraMeshPointInsertionProcess";
	}

	struct {
		std::unordered_map<std::size_t, std::size_t> nodes;
		std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
		std::unordered_map<std::size_t, std::list<std::size_t>> connect;
	} Graph;

	void TetrahedraMeshPointInsertionProcess::Execute() {

		Point insert_point;
		Graph mesh_cut;


		CalculateInsertPoint(mrModelPart, insert_point);
		CalculateCut(mrModelPart, insert_point, depth, mesh_cut);

	}

	void TetrahedraMeshPointInsertionProcess::CalculateInsertPoint(ModelPart & rModelPart, Element & rElement, Point<double,3> & rInsertPoint)
	{



		auto const& r_neighbours = rNode.GetValue(NEIGHBOUR_ELEMENTS);
		const std::size_t size = r_neighbours.size();
		rOptimumPoints.resize(size, ZeroVector(3));
		rWeights.resize(size);
		for (std::size_t i = 0; i < size; i++)
		{
			CalculateElementOptimumPosition(rNode, r_neighbours[i].GetGeometry(), rOptimumPoints[i]);
			auto quality = std::abs(r_neighbours[i].GetGeometry().Quality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH));
			if (quality > 1e-6)
				rWeights[i] = 1.00 / quality;
			else
				rWeights[i] = 1e6;
		}
	}

	void VisitElement(Element & rElem, std::list<Node<3>::Pointer> nodes_to_traverse, std::list<Element::Pointer> elems_to_traverse, Graph & mesh_cut, int travelHash) {
		int & this_hash = rElem->FastGetValue(TRAVEL_HASH);
		if(this_hash != travelHash) {
			this_hash = travelHash;
		}
	}

	void WalkGraph(ModelPart & rModelPart, Element & rRoot, std::list<Element::Pointer> elems_to_traverse, Graph & mesh_cut, int travelHash, std::size_t currentDepth, std::size_t maxDepth) {
		
		// If maximum depth has not been reached continue the walk, otherwise stop
		if(currentDepth < maxDepth) {
			std::list<Element::Pointer> elems_to_walk;

			auto root_id = root->Id();
			auto & neighbour_elements = NeighbourLocator();

			for(auto & neighbour: neighbour_elements) {
				auto & neighbour_hash = neighbour->GetValue(TRAVEL_HASH);
				auto neighbour_id = neighbour->Id();

				if(travel_hash != neighbour_hash) {
					neighbour_hash = travel_hash;
					mesh_cut.nodes.insert(std::make_pair<std::size_t, std::size_t>(neighbour_id, current_depth));
					mesh_cut.edges.insert(std::make_pair<std::pair<std::size_t, std::size_t>, std::size_t>(std::make_pair<std::size_t, std::size_t>(root_id, neighbour_id), 0));
					auto connect_list = mesh_cut.connect.find(root_id);
					if(connect_list != std::unordered_map::end) {
						mesh_cut.connect.insert(std::make_pair<std::size_t, std::list<std::size_t>>(root_id,std::list<std::size_t>));
						connect_list = mesh_cut.connect.find(root_id);
					}
					connect_list->second.append(neighbour_id);
					elems_to_walk.append(neighbour);
				}
			}
		}

		for(auto & neighbour: elems_to_walk) {
			WalkGraph(rModelPart, neighbour, mesh_cut, current_depth, depth);
		}
	}

	void TetrahedraMeshPointInsertionProcess::CalculateCut(ModelPart & rModelPart, Element & rElement, Point<double,3> & rInsertPoint, Graph & mesh_cut, std::size_t depth) {		
		std::list<Element::Pointer> elems_to_walk;

		auto & root_elements = FastPointLocator(rInsertPoint);
		int travel_hash = GenerateRandomHash();
		std::size_t current_depth = 0;

		for(auto & root: RootElements) {
			auto & root_hash = root->GetValue(TRAVEL_HASH);
			auto root_id = root->Id();
			
			// Add all possible roots not visited to the list of elements to traverse
			if(travel_hash != root_hash) {
				root_hash = travel_hash;
				mesh_cut.nodes.insert(std::make_pair<std::size_t, std::size_t>(root_id, current_depth));
				elems_to_walk.append(root);
			}
		}

		for(auto & root: elems_to_walk) {
			WalkGraph(rModelPart, root, mesh_cut, current_depth, depth);
		}
		
	}

	void TetrahedraMeshPointInsertionProcess::CalculateCut() {

	}


}  // namespace Kratos.
