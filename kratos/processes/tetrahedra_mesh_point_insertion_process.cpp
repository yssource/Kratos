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
#include <unordered_set>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/exception.h"
#include "includes/kratos_flags.h"
#include "processes/tetrahedra_mesh_point_insertion_process.h"
#include "processes/find_elements_neighbours_process.h"

namespace Kratos {

using Graph = TetrahedraMeshPointInsertionProcess::Graph;

TetrahedraMeshPointInsertionProcess::TetrahedraMeshPointInsertionProcess(ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart)
{
}

TetrahedraMeshPointInsertionProcess::~TetrahedraMeshPointInsertionProcess() 
{
}

/// Turn back information as a string.
std::string TetrahedraMeshPointInsertionProcess::Info() const
{
    return "TetrahedraMeshPointInsertionProcess";
}

void TetrahedraMeshPointInsertionProcess::Execute() 
{
    Point insert_point;
    Graph mesh_cut;

    Element::Pointer target_element;

    int depth = 6;

    FindElementalNeighboursProcess(mrModelPart, 3).Execute();

    CalculateInsertPoint(target_element, insert_point);
    CalculateCut(mrModelPart, target_element, insert_point, mesh_cut, depth);
}

/**
 * @brief Obtains the candidate coordinates to insert a point.
 * TODO: Now it only returns the center of the element, needs to be expanded.
 * 
 * @param rElement element in which to point will be inserted
 * @param rInsertPoint coordinates of the insertion
 */
void TetrahedraMeshPointInsertionProcess::CalculateInsertPoint(Element::Pointer & pElement, Point & rInsertPoint)
{
    rInsertPoint = pElement->GetGeometry().Center();
}


void VisitElement(Element & rElem, std::list<Node<3>::Pointer> nodesToTraverse, std::list<Element::Pointer> elemsToTraverse, Graph & meshCut, int travelHash)
{
    // int & this_hash = rElem->FastGetValue(TRAVEL_HASH);
    // if(this_hash != travelHash) {
    //     this_hash = travelHash;
    // }
}

void WalkGraph(ModelPart & rModelPart, Element::Pointer & pRootElem, Graph & meshCut, int travelHash, std::size_t currentDepth, std::size_t maxDepth)
{
    // // If maximum depth has not been reached continue the walk, otherwise stop
    // if(currentDepth < maxDepth) {
    //     std::list<Element::Pointer> elems_to_walk;

    //     auto root_id = root->Id();
    //     auto & neighbour_elements = NeighbourLocator();

    //     for(auto & neighbour: neighbour_elements) {
    //         auto & neighbour_hash = neighbour->GetValue(TRAVEL_HASH);
    //         auto neighbour_id = neighbour->Id();

    //         if(travel_hash != neighbour_hash) {
    //             neighbour_hash = travel_hash;
    //             mesh_cut.nodes.insert(std::make_pair<std::size_t, std::size_t>(neighbour_id, current_depth));
    //             mesh_cut.edges.insert(std::make_pair<std::pair<std::size_t, std::size_t>, std::size_t>(std::make_pair<std::size_t, std::size_t>(root_id, neighbour_id), 0));
    //             auto connect_list = mesh_cut.connect.find(root_id);
    //             if(connect_list != std::unordered_map::end) {
    //                 mesh_cut.connect.insert(std::make_pair<std::size_t, std::list<std::size_t>>(root_id,std::list<std::size_t>));
    //                 connect_list = mesh_cut.connect.find(root_id);
    //             }
    //             connect_list->second.append(neighbour_id);
    //             elems_to_walk.append(neighbour);
    //         }
    //     }
    // }

    // for(auto & neighbour: elems_to_walk) {
    //     WalkGraph(rModelPart, neighbour, mesh_cut, current_depth, depth);
    // }
}

/**
 * @brief To be implemented
 * 
 * @return int Random hash key
 */
int GenerateRandomHash() 
{
    return 4;
}

/**
 * @brief 
 * 
 * @param rModelPart 
 * @param rElement 
 * @param rInsertPoint 
 * @param meshCut 
 * @param depth 
 */
void TetrahedraMeshPointInsertionProcess::CalculateCut(ModelPart & rModelPart, Element::Pointer & pElement, Point & rInsertPoint, Graph & meshCut, std::size_t depth)
{		
    PointerVectorSet<Element> cut;
    PointerVectorSet<Element> roots;

    int travel_hash = GenerateRandomHash();
    
    std::size_t current_depth = 0;
    std::size_t max_depth = depth;

    int & root_hash = pElement->GetValue(STEP);
    std::size_t root_id = pElement->Id();
        
    // Add all possible roots not visited to the list of elements to traverse (usually 1)
    if(travel_hash != root_hash) {
        root_hash = travel_hash;
        cut.push_back(pElement);
        roots.push_back(pElement);
    }

    // Walk all the roots
    for(auto root = roots.ptr_begin(); root != roots.ptr_end(); root++) {
        WalkGraph(rModelPart, *root, meshCut, travel_hash, current_depth, max_depth);
    }	
}

} // namespace Kratos.
