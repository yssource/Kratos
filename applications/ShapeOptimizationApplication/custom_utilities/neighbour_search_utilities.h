// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef NEIGHBOUR_SEARCH_UTILITIES_H
#define NEIGHBOUR_SEARCH_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "shape_optimization_application.h"

// ==============================================================================

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
/** Detail class definition.

*/

class NeighbourSearchUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions
    // ==========================================================================
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;    
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;    

    /// Pointer definition of NeighbourSearchUtilities
    KRATOS_CLASS_POINTER_DEFINITION(NeighbourSearchUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NeighbourSearchUtilities( ModelPart& model_part )
    {
        // Create list of nodes for search tree
        NodeVector list_of_nodes(model_part.Nodes().size());
        int counter = 0;
        for(ModelPart::NodesContainerType::iterator node_it = model_part.NodesBegin(); node_it != model_part.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            list_of_nodes[counter++] = pnode;
        }

        // Create search tree
        mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(list_of_nodes.begin(), list_of_nodes.end(), mBucketSize));
    }

    /// Destructor.
    virtual ~NeighbourSearchUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    boost::python::dict SearchInRadius(NodeType& node_of_interest, double search_radius )
    {
        KRATOS_TRY;

        NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
        std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
        unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_of_interest,
                                                                         search_radius,
                                                                         neighbor_nodes.begin(),
                                                                         resulting_squared_distances.begin(),
                                                                         mMaxNumberOfNeighbors );

        boost::python::dict results;
        boost::python::list list_of_neighbours;
        boost::python::list list_of_neighbour_distances;

        for(int i=0; i<number_of_neighbors; i++)
            list_of_neighbours.append(neighbor_nodes[i]->Id());
        results["neighbour_node_ids"] = list_of_neighbours;

        for(int i=0; i<number_of_neighbors; i++)
            list_of_neighbour_distances.append( sqrt(resulting_squared_distances[i]) );
        results["neighbour_node_distances"] = list_of_neighbour_distances;

        return results;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    boost::python::dict SearchNearestPoint(NodeType& node_of_interest)
    {
        KRATOS_TRY;

        // double distance = 0.0;
        // Node<3>::Pointer nearest_point = mpSearchTree->SearchNearestPoint( node_of_interest, distance );

        boost::python::dict results;
        // results["nearest_node_id"] = nearest_point->Id();
        // results["nearest_node_distance"] = distance;

        return results;

        KRATOS_CATCH("");
    }

    // ==============================================================================

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
    virtual std::string Info() const
    {
        return "NeighbourSearchUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NeighbourSearchUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    unsigned int mBucketSize = 100;
    unsigned int mMaxNumberOfNeighbors = 100000;
    KDTree::Pointer mpSearchTree;    

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
//      NeighbourSearchUtilities& operator=(NeighbourSearchUtilities const& rOther);

    /// Copy constructor.
//      NeighbourSearchUtilities(NeighbourSearchUtilities const& rOther);


    ///@}

}; // Class NeighbourSearchUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // NEIGHBOUR_SEARCH_UTILITIES_H
