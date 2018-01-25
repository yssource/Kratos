// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Willer Matthias, https://github.com/matthiaswiller
//
// ==============================================================================

#ifndef RIGID_BODY_UTILITIES_H
#define RIGID_BODY_UTILITIES_H

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
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/normal_calculation_utils.h"
#include "spaces/ublas_space.h"
#include "shape_optimization_application.h"
#include "utilities/svd_utils.h"

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

class RigidBodyUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperVertexMorphing
    // KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidBodyUtilities( ModelPart& designSurface, Parameters optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size())
    {
        CreateListOfRigidNodes();
    }

    /// Destructor.
    virtual ~RigidBodyUtilities()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateListOfRigidNodes()
    {
        mListOfRigidNodes.resize(mNumberOfDesignVariables);
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            array_3d& nodal_variable = node_it.FastGetSolutionStepValue(rNodalVariable);
            if(nodal_variable[0] < -5.0)
            {
                NodeTypePointer pnode = *(node_it.base());
                mListOfRigidNodes[counter++] = pnode;
            }
        }
        mNumberOfRigidNodes = counter+1;
    }

    // --------------------------------------------------------------------------
    void FindRotationAndTranslationWithSVD()
    {
        // centroid A (member variable)
        Vector centroid_undeformed = ZeroVector(3);
        // centroid B (member variable)
        Vector centroid_deformed = ZeroVector(3);
        array_3d& coord;
        array_3d& def;
        for(auto& node_i : mRigidNodes.Nodes())
        {
            coord = node_i.Coordinates();
            def = node_i.FastGetSolutionStepValue( SHAPE_UPDATE );

            centroid_undeformed(0) += coord[0];
            centroid_undeformed(1) += coord[1];
            centroid_undeformed(2) += coord[2];

            centroid_deformed(0) += coord[0] + def[0];
            centroid_deformed(1) += coord[1] + def[1];
            centroid_deformed(2) += coord[2] + def[2];
        }
        centroid_undeformed(0) = centroid_undeformed(0)/mNumberOfRigidNodes;
        centroid_undeformed(1) = centroid_undeformed(1)/mNumberOfRigidNodes;
        centroid_undeformed(2) = centroid_undeformed(2)/mNumberOfRigidNodes;
        centroid_deformed(0) = centroid_deformed(0)/mNumberOfRigidNodes;
        centroid_deformed(1) = centroid_deformed(1)/mNumberOfRigidNodes;
        centroid_deformed(2) = centroid_deformed(2)/mNumberOfRigidNodes;

        // H matrix (3x3)
        Matrix H = ZeroMatrix(3,3);
        for (int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                for(auto& node_i : mListOfRigidNodes)
                {
                    array_3d& coord = node_i.Coordinates();
                    array_3d& def = node_i.FastGetSolutionStepValue( SHAPE_UPDATE );
                    H(i,j) += (coord[i] - centroid_undeformed(i))*(coord[j] + def[j] - centroid_deformed(j));
                }
            }
        }

        // U,V,S = svd(H)
        Matrix U = ZeroMatrix(3,3);
        Matrix V = ZeroMatrix(3,3);
        Matrix S = ZeroMatrix(3,3);

        SingularValueDecomposition( H&, U&, S&, V& );

        // R
        Matrix R;
        multiply(U,V.T,R&)

        if (determinant(R) < 0) // special reflection case
        {
            R(0,0) *= -1;
            R(0,1) *= -1;
            R(0,2) *= -1;
        }

        // t
        Vector t = ZeroVector(3);
        t = - R*centroid_undeformed + centroid_undeformed;

        print(R);
        print(t);
    }

    // --------------------------------------------------------------------------
    void CorrectDesignUpdateWithRigidBodyConstraints()
    {
        
    }

    // --------------------------------------------------------------------------
    // void InitializeMappingVariables()
    // {
    //     mMappingMatrix.resize(mNumberOfDesignVariables,mNumberOfDesignVariables);
    //     mMappingMatrix.clear();

    //     x_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
    //     y_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
    //     z_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);

    //     x_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
    //     y_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
    //     z_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
    // }

    // // --------------------------------------------------------------------------
    // void AssignMappingIds()
    // {
    //     unsigned int i = 0;
    //     for(auto& node_i : mrDesignSurface.Nodes())
    //         node_i.SetValue(MAPPING_ID,i++);
    // }

    // // --------------------------------------------------------------------------
    // void ComputeMappingMatrix()
    // {
    //     boost::timer timer;
    //     std::cout << "> Computing mapping matrix to perform mapping..." << std::endl;

    //     CreateSearchTreeWithAllNodesOnDesignSurface();
    //     ComputeEntriesOfMappingMatrix();

    //     std::cout << "> Mapping matrix computed in: " << timer.elapsed() << " s" << std::endl;
    // }

    // // --------------------------------------------------------------------------
    // void CreateSearchTreeWithAllNodesOnDesignSurface()
    // {
    //     mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mListOfNodesOfDesignSurface.begin(), mListOfNodesOfDesignSurface.end(), mBucketSize));
    // }

    // // --------------------------------------------------------------------------
    // void ComputeEntriesOfMappingMatrix()
    // {
    //     for(auto& node_i : mrDesignSurface.Nodes())
    //     {
    //         NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
    //         std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
    //         unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
    //                                                                          mFilterRadius,
    //                                                                          neighbor_nodes.begin(),
    //                                                                          resulting_squared_distances.begin(),
    //                                                                          mMaxNumberOfNeighbors );

    //         std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
    //         double sum_of_weights = 0.0;

    //         ThrowWarningIfMaxNodeNeighborsReached( node_i, number_of_neighbors );
    //         ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
    //         FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void ThrowWarningIfMaxNodeNeighborsReached( ModelPart::NodeType& given_node, unsigned int number_of_neighbors )
    // {
    //     if(number_of_neighbors >= mMaxNumberOfNeighbors)
    //         std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << mMaxNumberOfNeighbors << " nodes) reached!" << std::endl;
    // }

    // // --------------------------------------------------------------------------
    // virtual void ComputeWeightForAllNeighbors(  ModelPart::NodeType& design_node,
    //                                     NodeVector& neighbor_nodes,
    //                                     unsigned int number_of_neighbors,
    //                                     std::vector<double>& list_of_weights,
    //                                     double& sum_of_weights )
    // {
    //     for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
    //     {
    //         ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
    //         double weight = 0;

    //         list_of_weights[neighbor_itr] = weight;
    //         sum_of_weights += weight;
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void FillMappingMatrixWithWeights(  ModelPart::NodeType& design_node,
    //                                     NodeVector& neighbor_nodes,
    //                                     unsigned int number_of_neighbors,
    //                                     std::vector<double>& list_of_weights,
    //                                     double& sum_of_weights )
    // {
    //     unsigned int row_id = design_node.GetValue(MAPPING_ID);
    //     for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
    //     {
    //         ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
    //         int collumn_id = neighbor_node.GetValue(MAPPING_ID);

    //         double weight = list_of_weights[neighbor_itr] / sum_of_weights;
    //         mMappingMatrix.push_back(row_id,collumn_id,weight);
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    // {
    //     boost::timer mapping_time;
    //     std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

    //     RecomputeMappingMatrixIfGeometryHasChanged();
    //     PrepareVectorsForMappingToDesignSpace( rNodalVariable );
    //     if (mConsistentBackwardMapping)
    //         MultiplyVectorsWithConsistentBackwardMappingMatrix();
    //     else
    //         MultiplyVectorsWithTransposeMappingMatrix();
    //     AssignResultingDesignVectorsToNodalVariable( rNodalVariableInDesignSpace );

    //     std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    // }

    // // --------------------------------------------------------------------------
    // void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    // {
    //     boost::timer mapping_time;
    //     std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

    //     RecomputeMappingMatrixIfGeometryHasChanged();
    //     PrepareVectorsForMappingToGeometrySpace( rNodalVariable );
    //     MultiplyVectorsWithMappingMatrix();
    //     AssignResultingGeometryVectorsToNodalVariable( rNodalVariableInGeometrySpace );

    //     std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    // }

    // // --------------------------------------------------------------------------
    // void RecomputeMappingMatrixIfGeometryHasChanged()
    // {
    //     if(HasGeometryChanged())
    //     {
    //         InitializeComputationOfMappingMatrix();
    //         ComputeMappingMatrix();
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void PrepareVectorsForMappingToDesignSpace( const Variable<array_3d> &rNodalVariable )
    // {
    //     x_variables_in_design_space.clear();
    //     y_variables_in_design_space.clear();
    //     z_variables_in_design_space.clear();
    //     x_variables_in_geometry_space.clear();
    //     y_variables_in_geometry_space.clear();
    //     z_variables_in_geometry_space.clear();

    //     for(auto& node_i : mrDesignSurface.Nodes())
    //     {
    //         int i = node_i.GetValue(MAPPING_ID);
    //         array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
    //         x_variables_in_geometry_space[i] = nodal_variable[0];
    //         y_variables_in_geometry_space[i] = nodal_variable[1];
    //         z_variables_in_geometry_space[i] = nodal_variable[2];
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void PrepareVectorsForMappingToGeometrySpace( const Variable<array_3d> &rNodalVariable )
    // {
    //     x_variables_in_design_space.clear();
    //     y_variables_in_design_space.clear();
    //     z_variables_in_design_space.clear();
    //     x_variables_in_geometry_space.clear();
    //     y_variables_in_geometry_space.clear();
    //     z_variables_in_geometry_space.clear();

    //     for(auto& node_i : mrDesignSurface.Nodes())
    //     {
    //         int i = node_i.GetValue(MAPPING_ID);
    //         array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
    //         x_variables_in_design_space[i] = nodal_variable[0];
    //         y_variables_in_design_space[i] = nodal_variable[1];
    //         z_variables_in_design_space[i] = nodal_variable[2];
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void MultiplyVectorsWithTransposeMappingMatrix()
    // {
    //     SparseSpaceType::TransposeMult(mMappingMatrix,x_variables_in_geometry_space,x_variables_in_design_space);
    //     SparseSpaceType::TransposeMult(mMappingMatrix,y_variables_in_geometry_space,y_variables_in_design_space);
    //     SparseSpaceType::TransposeMult(mMappingMatrix,z_variables_in_geometry_space,z_variables_in_design_space);
    // }

    // // --------------------------------------------------------------------------
    // void MultiplyVectorsWithConsistentBackwardMappingMatrix()
    // {
    //     // for the case of matching grids in geometry and design space, use the forward mapping matrix
    //     noalias(x_variables_in_design_space) = prod(mMappingMatrix,x_variables_in_geometry_space);
    //     noalias(y_variables_in_design_space) = prod(mMappingMatrix,y_variables_in_geometry_space);
    //     noalias(z_variables_in_design_space) = prod(mMappingMatrix,z_variables_in_geometry_space);
    // }

    // // --------------------------------------------------------------------------
    // void MultiplyVectorsWithMappingMatrix()
    // {
    //     noalias(x_variables_in_geometry_space) = prod(mMappingMatrix,x_variables_in_design_space);
    //     noalias(y_variables_in_geometry_space) = prod(mMappingMatrix,y_variables_in_design_space);
    //     noalias(z_variables_in_geometry_space) = prod(mMappingMatrix,z_variables_in_design_space);
    // }

    // // --------------------------------------------------------------------------
    // void AssignResultingDesignVectorsToNodalVariable( const Variable<array_3d> &rNodalVariable )
    // {
    //     for(auto& node_i : mrDesignSurface.Nodes())
    //     {
    //         int i = node_i.GetValue(MAPPING_ID);

    //         Vector node_vector = ZeroVector(3);
    //         node_vector(0) = x_variables_in_design_space[i];
    //         node_vector(1) = y_variables_in_design_space[i];
    //         node_vector(2) = z_variables_in_design_space[i];
    //         node_i.FastGetSolutionStepValue(rNodalVariable) = node_vector;
    //     }
    // }

    // // --------------------------------------------------------------------------
    // void AssignResultingGeometryVectorsToNodalVariable( const Variable<array_3d> &rNodalVariable )
    // {
    //     for(auto& node_i : mrDesignSurface.Nodes())
    //     {
    //         int i = node_i.GetValue(MAPPING_ID);

    //         Vector node_vector = ZeroVector(3);
    //         node_vector(0) = x_variables_in_geometry_space[i];
    //         node_vector(1) = y_variables_in_geometry_space[i];
    //         node_vector(2) = z_variables_in_geometry_space[i];
    //         node_i.FastGetSolutionStepValue(rNodalVariable) = node_vector;
    //     }
    // }

    // // --------------------------------------------------------------------------
    // bool HasGeometryChanged()
    // {
    //     double sumOfAllCoordinates = 0.0;
    //     for(auto& node_i : mrDesignSurface.Nodes())
    //     {
    //         array_3d& coord = node_i.Coordinates();
    //         sumOfAllCoordinates += coord[0] + coord[1] + coord[2];
    //     }

    //     if (mControlSum == sumOfAllCoordinates)
    //         return false;
    //     else
    //     {
    //         mControlSum = sumOfAllCoordinates;
    //         return true;
    //     }
    // }

    // // --------------------------------------------------------------------------
    // virtual void InitializeComputationOfMappingMatrix()
    // {
    //     mpSearchTree.reset();
    //     mMappingMatrix.clear();
    // }

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
        return "RigidBodyUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RigidBodyUtilities";
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

    // // ==============================================================================
    // // Initialized by class constructor
    // // ==============================================================================
    ModelPart& mrDesignSurface;

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

    // // ==============================================================================
    // // Initialized by class constructor
    // // ==============================================================================
    const unsigned int mNumberOfDesignVariables;
    unsigned int mNumberOfRigidNodes;
    std::string mFilterType;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;
    bool mConsistentBackwardMapping;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesOfDesignSurface;
    NodeVector mListOfRigidNodes;
    KDTree::Pointer mpSearchTree;

    // ==============================================================================
    // Variables for mapping
    // ==============================================================================
    SparseMatrixType mMappingMatrix;
    Vector x_variables_in_design_space, y_variables_in_design_space, z_variables_in_design_space;
    Vector x_variables_in_geometry_space, y_variables_in_geometry_space, z_variables_in_geometry_space;
    double mControlSum = 0.0;

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
//      RigidBodyUtilities& operator=(RigidBodyUtilities const& rOther);

    /// Copy constructor.
//      RigidBodyUtilities(RigidBodyUtilities const& rOther);


    ///@}

}; // Class RigidBodyUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // RIGID_BODY_UTILITIES_H
