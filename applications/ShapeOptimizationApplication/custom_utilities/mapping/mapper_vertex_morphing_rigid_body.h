// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_RIGID_BODY_H
#define MAPPER_VERTEX_MORPHING_RIGID_BODY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

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
#include "filter_function.h"
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

class MapperVertexMorphingRigidBody
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
    typedef UblasSpace<double, CompressedMatrix, Vector> CompressedSpaceType;
    typedef CompressedSpaceType::MatrixType CompressedMatrixType;
    typedef CompressedSpaceType::VectorType VectorType;

    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<CompressedSpaceType, DenseSpaceType > CompressedLinearSolverType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperVertexMorphingRigidBody
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingRigidBody);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingRigidBody( ModelPart& designSurface, Parameters mapper_settings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size()),
          mFilterType( mapper_settings["filter_function_type"].GetString() ),
          mFilterRadius( mapper_settings["filter_radius"].GetDouble() ),
          mMaxNumberOfNeighbors( mapper_settings["max_nodes_in_filter_radius"].GetInt() ),
          mConsistentBackwardMapping (mapper_settings["consistent_mapping_to_geometry_space"].GetBool() ),
          mMapperSettings( mapper_settings )
    {
        CreateListOfNodesOfDesignSurface();
        CreateFilterFunction();
        InitializeMappingVariables();
        AssignMappingIds();
    }

    /// Destructor.
    virtual ~MapperVertexMorphingRigidBody()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateListOfNodesOfDesignSurface()
    {
        mListOfNodesOfDesignSurface.resize(mNumberOfDesignVariables);
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesOfDesignSurface[counter++] = pnode;
        }
    }

    // --------------------------------------------------------------------------
    void CreateFilterFunction()
    {
        mpFilterFunction = boost::shared_ptr<FilterFunction>(new FilterFunction(mFilterType, mFilterRadius));
    }

    // --------------------------------------------------------------------------
    void InitializeMappingVariables()
    {
        mMappingMatrix.resize(mNumberOfDesignVariables,mNumberOfDesignVariables,false);
        mMappingMatrix.clear();

        x_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
        y_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
        z_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);

        x_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
        y_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
        z_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
    }

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for(auto& node_i : mrDesignSurface.Nodes())
            node_i.SetValue(MAPPING_ID,i++);
    }

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        boost::timer timer;
        std::cout << "> Computing mapping matrix to perform mapping..." << std::endl;

        CreateSearchTreeWithAllNodesOnDesignSurface();
        ComputeEntriesOfMappingMatrix();

        std::cout << "> Mapping matrix computed in: " << timer.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesOnDesignSurface()
    {
        mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mListOfNodesOfDesignSurface.begin(), mListOfNodesOfDesignSurface.end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeEntriesOfMappingMatrix()
    {
        for(auto& node_i : mrDesignSurface.Nodes())
        {
            NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
            std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                             mFilterRadius,
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(),
                                                                             mMaxNumberOfNeighbors );

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;

            ThrowWarningIfMaxNodeNeighborsReached( node_i, number_of_neighbors );
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
        }
    }

    // --------------------------------------------------------------------------
    void ThrowWarningIfMaxNodeNeighborsReached( ModelPart::NodeType& given_node, unsigned int number_of_neighbors )
    {
        if(number_of_neighbors >= mMaxNumberOfNeighbors)
            std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << mMaxNumberOfNeighbors << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    virtual void ComputeWeightForAllNeighbors(  ModelPart::NodeType& design_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double weight = mpFilterFunction->compute_weight( design_node.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
        }
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& design_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        unsigned int row_id = design_node.GetValue(MAPPING_ID);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            int collumn_id = neighbor_node.GetValue(MAPPING_ID);

            double weight = list_of_weights[neighbor_itr] / sum_of_weights;
            mMappingMatrix.insert_element(row_id,collumn_id,weight);
        }
    }

    // --------------------------------------------------------------------------
    void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

        RecomputeMappingMatrixIfGeometryHasChanged();
        PrepareVectorsForMappingToDesignSpace( rNodalVariable );
        if (mConsistentBackwardMapping)
            MultiplyVectorsWithConsistentBackwardMappingMatrix();
        else
            MultiplyVectorsWithTransposeMappingMatrix();
        AssignResultingDesignVectorsToNodalVariable( rNodalVariableInDesignSpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // // --------------------------------------------------------------------------
    // void MapToDesignSpaceWithRigidCorrection( const Variable<array_3d> &rNodalVariable,
    //                                           const Variable<array_3d> &rNodalVariableInDesignSpace,
    //                                           CompressedLinearSolverType::Pointer linear_solver )
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

    //     CorrectDesignUpdateWithRigidBodyConstraints( linear_solver );
    //     AssignResultingDesignVectorsToNodalVariable( rNodalVariableInDesignSpace );

    //     std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    // }

    // --------------------------------------------------------------------------
    void MapToGeometrySpaceWithRigidCorrection( const Variable<array_3d> &rNodalVariable,
                                                const Variable<array_3d> &rNodalVariableInGeometrySpace,
                                                CompressedLinearSolverType::Pointer linear_solver,
                                                boost::python::list rigidRegions,
                                                boost::python::list fixedRegions )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

        mrFixedRegions = fixedRegions;
        mrRigidRegions = rigidRegions;
        CreateListOfRigidNodes();
        CreateListOfFixedNodes();
        CreateListOfActiveNodes();

        RecomputeMappingMatrixIfGeometryHasChanged();
        PrepareVectorsForMappingToGeometrySpace( rNodalVariable );
        MultiplyVectorsWithMappingMatrix();
        AssignResultingGeometryVectorsToNodalVariable( rNodalVariableInGeometrySpace );

        // CorrectDesignUpdateWithRigidBodyConstraints( linear_solver );
        CorrectDesignUpdateWithRigidBodyConstraintsFast( linear_solver );
        MultiplyVectorsWithMappingMatrix();
        AssignResultingGeometryVectorsToNodalVariable( rNodalVariableInGeometrySpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

        RecomputeMappingMatrixIfGeometryHasChanged();
        PrepareVectorsForMappingToGeometrySpace( rNodalVariable );
        MultiplyVectorsWithMappingMatrix();
        AssignResultingGeometryVectorsToNodalVariable( rNodalVariableInGeometrySpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void RecomputeMappingMatrixIfGeometryHasChanged()
    {
        if(HasGeometryChanged())
        {
            InitializeComputationOfMappingMatrix();
            ComputeMappingMatrix();
        }
    }

    // --------------------------------------------------------------------------
    void PrepareVectorsForMappingToDesignSpace( const Variable<array_3d> &rNodalVariable )
    {
        x_variables_in_design_space.clear();
        y_variables_in_design_space.clear();
        z_variables_in_design_space.clear();
        x_variables_in_geometry_space.clear();
        y_variables_in_geometry_space.clear();
        z_variables_in_geometry_space.clear();

        for(auto& node_i : mrDesignSurface.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
            x_variables_in_geometry_space[i] = nodal_variable[0];
            y_variables_in_geometry_space[i] = nodal_variable[1];
            z_variables_in_geometry_space[i] = nodal_variable[2];
        }
    }

    // --------------------------------------------------------------------------
    void PrepareVectorsForMappingToGeometrySpace( const Variable<array_3d> &rNodalVariable )
    {
        x_variables_in_design_space.clear();
        y_variables_in_design_space.clear();
        z_variables_in_design_space.clear();
        x_variables_in_geometry_space.clear();
        y_variables_in_geometry_space.clear();
        z_variables_in_geometry_space.clear();

        for(auto& node_i : mrDesignSurface.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
            x_variables_in_design_space[i] = nodal_variable[0];
            y_variables_in_design_space[i] = nodal_variable[1];
            z_variables_in_design_space[i] = nodal_variable[2];
        }
    }

    // --------------------------------------------------------------------------
    void MultiplyVectorsWithTransposeMappingMatrix()
    {
        CompressedSpaceType::TransposeMult(mMappingMatrix,x_variables_in_geometry_space,x_variables_in_design_space);
        CompressedSpaceType::TransposeMult(mMappingMatrix,y_variables_in_geometry_space,y_variables_in_design_space);
        CompressedSpaceType::TransposeMult(mMappingMatrix,z_variables_in_geometry_space,z_variables_in_design_space);
    }

    // --------------------------------------------------------------------------
    void MultiplyVectorsWithConsistentBackwardMappingMatrix()
    {
        // for the case of matching grids in geometry and design space, use the forward mapping matrix
        noalias(x_variables_in_design_space) = prod(mMappingMatrix,x_variables_in_geometry_space);
        noalias(y_variables_in_design_space) = prod(mMappingMatrix,y_variables_in_geometry_space);
        noalias(z_variables_in_design_space) = prod(mMappingMatrix,z_variables_in_geometry_space);
    }

    // --------------------------------------------------------------------------
    void MultiplyVectorsWithMappingMatrix()
    {
        noalias(x_variables_in_geometry_space) = prod(mMappingMatrix,x_variables_in_design_space);
        noalias(y_variables_in_geometry_space) = prod(mMappingMatrix,y_variables_in_design_space);
        noalias(z_variables_in_geometry_space) = prod(mMappingMatrix,z_variables_in_design_space);
    }

    // --------------------------------------------------------------------------
    void AssignResultingDesignVectorsToNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for(auto& node_i : mrDesignSurface.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            Vector node_vector = ZeroVector(3);
            node_vector(0) = x_variables_in_design_space[i];
            node_vector(1) = y_variables_in_design_space[i];
            node_vector(2) = z_variables_in_design_space[i];
            node_i.FastGetSolutionStepValue(rNodalVariable) = node_vector;
        }
    }

    // --------------------------------------------------------------------------
    void AssignResultingGeometryVectorsToNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for(auto& node_i : mrDesignSurface.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            Vector node_vector = ZeroVector(3);
            node_vector(0) = x_variables_in_geometry_space[i];
            node_vector(1) = y_variables_in_geometry_space[i];
            node_vector(2) = z_variables_in_geometry_space[i];
            node_i.FastGetSolutionStepValue(rNodalVariable) = node_vector;
        }
    }

    // --------------------------------------------------------------------------
    bool HasGeometryChanged()
    {
        double sumOfAllCoordinates = 0.0;
        for(auto& node_i : mrDesignSurface.Nodes())
        {
            array_3d& coord = node_i.Coordinates();
            sumOfAllCoordinates += coord[0] + coord[1] + coord[2];
        }

        if (mControlSum == sumOfAllCoordinates)
            return false;
        else
        {
            mControlSum = sumOfAllCoordinates;
            return true;
        }
    }

    // --------------------------------------------------------------------------
    virtual void InitializeComputationOfMappingMatrix()
    {
        mpSearchTree.reset();
        mMappingMatrix.clear();
    }


    // --------------------------------------------------------------------------
    void CreateListOfFixedNodes()
    {
        mListOfFixedNodes.clear();
        // Loop over all regions that have to be fixed
        for (unsigned int regionNumber = 0; regionNumber < len(mrFixedRegions); regionNumber++)
        {
            ModelPart& fixed_region = extractModelPart( mrFixedRegions[regionNumber] );
            for (ModelPart::NodesContainerType::iterator node_it = fixed_region.NodesBegin(); node_it != fixed_region.NodesEnd(); ++node_it)
            {
                NodeTypePointer pnode = *(node_it.base());
                mListOfFixedNodes.push_back( pnode );
            }

            // mFix_in_X = mMapperSettings[regionNumber]["fix_in_X"].GetBool();
            // mFix_in_Y = mMapperSettings[regionNumber]["fix_in_Y"].GetBool();
            // mFix_in_Z = mMapperSettings[regionNumber]["fix_in_Z"].GetBool();
            mFix_in_X = true;
            mFix_in_Y = true;
            mFix_in_Z = true;
        }

        // writeVectorToFile(mListOfFixedNodes, "listOfFixedNodes.txt");
    }


    // --------------------------------------------------------------------------
    void CreateListOfRigidNodes()
    {
        mListOfRigidNodes.clear();
        // Loop over all regions that a rigid motion has to be enforced on
        for (unsigned int regionNumber = 0; regionNumber < len(mrRigidRegions); regionNumber++)
        {
            ModelPart& rigid_region = extractModelPart( mrRigidRegions[regionNumber] );
            for (ModelPart::NodesContainerType::iterator node_it = rigid_region.NodesBegin(); node_it != rigid_region.NodesEnd(); ++node_it)
            {
                NodeTypePointer pnode = *(node_it.base());
                mListOfRigidNodes.push_back( pnode );
            }
        }
        mNumberOfRigidNodes = mListOfRigidNodes.size();

        // writeVectorToFile(mListOfRigidNodes, "listOfRigidNodes.txt");
    }

    // --------------------------------------------------------------------------
    void CreateListOfActiveNodes()
    {
        double important_radius = mFilterRadius;

        std::vector<bool> boolVectorOfActiveNodes;
        boolVectorOfActiveNodes.resize(mNumberOfDesignVariables,false);

        // Loop over all rigid regions
        for (unsigned int regionNumber = 0; regionNumber < len(mrRigidRegions); regionNumber++)
        {
            ModelPart& rigidRegion = extractModelPart( mrRigidRegions[regionNumber] );

            // Loop over all nodes in specified sub-model part
            for(auto& node_i : rigidRegion.Nodes())
            {
                NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
                std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
                unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                                important_radius,
                                                                                neighbor_nodes.begin(),
                                                                                resulting_squared_distances.begin(),
                                                                                mMaxNumberOfNeighbors );

                ThrowWarningIfMaxNodeNeighborsReached( node_i, number_of_neighbors );

                for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
                {
                    ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                    int neighbor_id = neighbor_node.GetValue( MAPPING_ID );
                    boolVectorOfActiveNodes[neighbor_id] = true;
                }
            }
        }

        // Loop over all fixed regions
        for (unsigned int regionNumber = 0; regionNumber < len(mrFixedRegions); regionNumber++)
        {
            ModelPart& fixedRegion = extractModelPart( mrFixedRegions[regionNumber] );

            // Loop over all nodes in specified sub-model part
            for(auto& node_i : fixedRegion.Nodes())
            {
                NodeVector neighbor_nodes( mMaxNumberOfNeighbors );
                std::vector<double> resulting_squared_distances( mMaxNumberOfNeighbors );
                unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                                important_radius,
                                                                                neighbor_nodes.begin(),
                                                                                resulting_squared_distances.begin(),
                                                                                mMaxNumberOfNeighbors );

                ThrowWarningIfMaxNodeNeighborsReached( node_i, number_of_neighbors );

                for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
                {
                    ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                    int neighbor_id = neighbor_node.GetValue( MAPPING_ID );
                    boolVectorOfActiveNodes[neighbor_id] = true;
                }
            }
        }


        // create lists of important nodes
        mListOfActiveNodes.clear();
        for( auto& node_i : mListOfNodesOfDesignSurface )
        {
            int id = node_i->GetValue( MAPPING_ID );

            if ( boolVectorOfActiveNodes[id] )
                mListOfActiveNodes.push_back( node_i );
        }

        // writeVectorToFile(mListOfActiveNodes, "listOfActiveNodes.txt");
    }

    // --------------------------------------------------------------------------
    void computeTranslationVectorAndRotationMatrix()
    {
        //initialization
        mRotationMatrix = ZeroMatrix(3,3);
        mTranslationVector = ZeroVector(3);

        if(mListOfRigidNodes.size() > 0)
        {
            boost::timer svd_time;
            std::cout << "\n> Starting to compute SVD ..." << std::endl;

            Vector centroid_undeformed = ZeroVector(3); // centroid A
            Vector centroid_deformed = ZeroVector(3); // centroid B

            // compute centroid deformed/undeformed
            for( auto& node_i : mListOfRigidNodes )
            {
                Vector coord = node_i->Coordinates();
                Vector def = node_i->FastGetSolutionStepValue( SHAPE_UPDATE );

                centroid_undeformed += coord;
                centroid_deformed += coord + def;
            }
            centroid_undeformed = centroid_undeformed/mNumberOfRigidNodes;
            centroid_deformed = centroid_deformed/mNumberOfRigidNodes;

            // H matrix (3x3)
            Matrix H = ZeroMatrix(3,3);

            // assemble H-matrix
            for( auto& node_i : mListOfRigidNodes )
            {
                Vector coord = node_i->Coordinates();
                Vector def = node_i->FastGetSolutionStepValue( SHAPE_UPDATE );

                H += outer_prod( (coord - centroid_undeformed) , (coord + def - centroid_deformed) );
            }

            // U,V,S = svd(H)
            Matrix U = ZeroMatrix(3,3);
            Matrix V = ZeroMatrix(3,3);
            Matrix S = ZeroMatrix(3,3);

            SVDUtils<double>::JacobiSingularValueDecomposition(H, U, S, V);

            // mRotationMatrix [R], note that V is computed as transpose in the above svd algorithm
            mRotationMatrix = prod(trans(V),trans(U));


            // for ( int i_itr = 0; i_itr<3; i_itr++)
            //     if(std::abs(mRotationMatrix(i_itr,i_itr))<1e-10)
            //         mRotationMatrix(i_itr,i_itr) = 1.0;
            // for ( int i_itr = 0; i_itr<3; i_itr++)
            //     for ( int j_itr = 0; j_itr<3; j_itr++)
            //         if(std::abs(mRotationMatrix(i_itr,j_itr))<1e-10)

            double detR = MathUtils<double>::Det(mRotationMatrix);

            if ( detR < 0 ) // special reflection case
            {
                // Remove this warning, when the method works flawlessly!
                std::cout << "> Warning: det(R) < 0. Entering special reflection case!" << std::endl;

                mRotationMatrix(0,2) *= -1;
                mRotationMatrix(1,2) *= -1;
                mRotationMatrix(2,2) *= -1;
            }

            // mTranslationVector [t]
            mTranslationVector = - prod(mRotationMatrix, centroid_undeformed) + centroid_deformed;

            std::cout << "> Time needed for SVD: " << svd_time.elapsed() << " s" << std::endl;
        }
    }

    // --------------------------------------------------------------------------
    void CorrectDesignUpdateWithRigidBodyConstraints( CompressedLinearSolverType::Pointer linear_solver )
    {
        computeTranslationVectorAndRotationMatrix();

        boost::timer assembling_time;
        std::cout << "\n> Starting to assemble modifiedMappingMatrix..." << std::endl;

        // compute modified mapping matrix
        int numberOfRowsInModifiedSystem = mNumberOfDesignVariables + mListOfRigidNodes.size() + mListOfFixedNodes.size();
        CompressedMatrixType modifiedMappingMatrix( numberOfRowsInModifiedSystem , mNumberOfDesignVariables );

        Vector x_variables_modified, y_variables_modified, z_variables_modified;
        x_variables_modified.resize( numberOfRowsInModifiedSystem , 0.0 );
        y_variables_modified.resize( numberOfRowsInModifiedSystem , 0.0 );
        z_variables_modified.resize( numberOfRowsInModifiedSystem , 0.0 );

        double penalty_factor = 100000;

        // Using identity matrix for initial ds
        int counter = 0;
        for( auto & node_i : mrDesignSurface.Nodes() )
        {
            modifiedMappingMatrix.insert_element(counter, counter, 1.0);

            // modified vectors
            x_variables_modified[counter] = x_variables_in_design_space[counter];
            y_variables_modified[counter] = y_variables_in_design_space[counter];
            z_variables_modified[counter] = z_variables_in_design_space[counter];

            counter ++;
        }

        // rigid body transformation
        for( auto& node_i : mListOfRigidNodes )
        {
            // Get node information
            int i = node_i->GetValue( MAPPING_ID );

            // Copy respective row from mMappingMatrix to modifiedMappingMatrix
            row( modifiedMappingMatrix , counter ) = penalty_factor * row( mMappingMatrix , i );

            // compute rigid body movement of rigid nodes and store it in modified vector
            Vector coord = node_i->Coordinates();

            Vector modifiedDeformation = prod(mRotationMatrix,coord) + mTranslationVector - coord;

            x_variables_modified[counter] = penalty_factor*modifiedDeformation(0);
            y_variables_modified[counter] = penalty_factor*modifiedDeformation(1);
            z_variables_modified[counter] = penalty_factor*modifiedDeformation(2);

            counter ++;
        }

        // fixed points
        for( auto& node_i : mListOfFixedNodes )
        {
            // Get node information
            int i = node_i->GetValue( MAPPING_ID );

            // Copy respective row from mMappingMatrix to modifiedMappingMatrix
            row( modifiedMappingMatrix , counter ) = penalty_factor * row( mMappingMatrix , i );

            x_variables_modified[counter] = 0;
            y_variables_modified[counter] = 0;
            z_variables_modified[counter] = 0;

            if (!mFix_in_X)
                x_variables_modified[counter] = penalty_factor * x_variables_in_geometry_space[i];
            if (!mFix_in_Y)
                y_variables_modified[counter] = penalty_factor * y_variables_in_geometry_space[i];
            if (!mFix_in_Z)
                z_variables_modified[counter] = penalty_factor * z_variables_in_geometry_space[i];

            counter ++;
        }

        std::cout << "> Time needed for assembling modifiedMappingMatrix: " << assembling_time.elapsed() << " s" << std::endl;

        // write mapping matrix

        // writeMatrixToFile( modifiedMappingMatrix,
        //                    numberOfRowsInModifiedSystem,
        //                    mNumberOfDesignVariables,
        //                    "modifiedMappingMatrixFull.txt");

        // writeVectorToFile ( x_variables_modified,
        //                     "x_variables_modified_full.txt");
        // writeVectorToFile ( y_variables_modified,
        //                     "y_variables_modified_full.txt");
        // writeVectorToFile ( z_variables_modified,
        //                     "z_variables_modified_full.txt");

        // writeMatrixToFile( mMappingMatrix,
        //                    mNumberOfDesignVariables,
        //                    mNumberOfDesignVariables,
        //                    "originalMappingMatrix.txt");

        boost::timer solving_time;
        std::cout << "\n> Starting to solve modified System..." << std::endl;

        x_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
        y_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
        z_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);

        linear_solver->Solve(modifiedMappingMatrix, x_variables_in_design_space, x_variables_modified);
        linear_solver->Solve(modifiedMappingMatrix, y_variables_in_design_space, y_variables_modified);
        linear_solver->Solve(modifiedMappingMatrix, z_variables_in_design_space, z_variables_modified);

        // writeVectorToFile ( x_variables_in_design_space,
        //                     "x_variables_in_design_space_full.txt");
        // writeVectorToFile ( y_variables_in_design_space,
        //                     "y_variables_in_design_space_full.txt");
        // writeVectorToFile ( z_variables_in_design_space,
        //                     "z_variables_in_design_space_full.txt");

        std::cout << "> Time needed for solving modified System: " << solving_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    //
    // --------------------------------------------------------------------------
    void CorrectDesignUpdateWithRigidBodyConstraintsFast( CompressedLinearSolverType::Pointer linear_solver )
    {
        computeTranslationVectorAndRotationMatrix();

        boost::timer assembling_time;
        std::cout << "\n> Starting to assemble modifiedMappingMatrix..." << std::endl;

        // compute modified mapping matrix
        unsigned int numberOfRowsInModifiedSystem = mListOfActiveNodes.size() + mListOfRigidNodes.size() + mListOfFixedNodes.size();
        unsigned int numberOfColsInModifiedSystem = mListOfActiveNodes.size();

        CompressedMatrixType modifiedMappingMatrix( numberOfRowsInModifiedSystem , numberOfColsInModifiedSystem );

        Vector x_variables_modified, y_variables_modified, z_variables_modified;
        x_variables_modified.resize( numberOfRowsInModifiedSystem , 0.0 );
        y_variables_modified.resize( numberOfRowsInModifiedSystem , 0.0 );
        z_variables_modified.resize( numberOfRowsInModifiedSystem , 0.0 );

        double penalty_factor = 100000;

        // Using identity matrix for initial ds
        int counter_row = 0;
        for( auto& node_i : mListOfActiveNodes )
        {
            int i = node_i->GetValue( MAPPING_ID );
            modifiedMappingMatrix.insert_element(counter_row, counter_row, 1.0);

            // modified vectors
            x_variables_modified[counter_row] = x_variables_in_design_space[i];
            y_variables_modified[counter_row] = y_variables_in_design_space[i];
            z_variables_modified[counter_row] = z_variables_in_design_space[i];

            counter_row ++;
        }

        // rigid body transformation
        for( auto& node_i : mListOfRigidNodes )
        {
            // Get node information
            int i = node_i->GetValue( MAPPING_ID );
            int counter_col = 0;
            for( auto& node_j : mListOfActiveNodes )
            {
                int j = node_j->GetValue( MAPPING_ID );
                double val = mMappingMatrix(i,j);
                if (val != 0.0)
                    modifiedMappingMatrix.insert_element( counter_row, counter_col, penalty_factor*val );
                counter_col ++;
            }

            // compute rigid body movement of rigid nodes and store it in modified vector
            Vector coord = node_i->Coordinates();
            Vector modifiedDeformation = prod(mRotationMatrix,coord) + mTranslationVector - coord;

            x_variables_modified[counter_row] = penalty_factor*modifiedDeformation(0);
            y_variables_modified[counter_row] = penalty_factor*modifiedDeformation(1);
            z_variables_modified[counter_row] = penalty_factor*modifiedDeformation(2);
            counter_row ++;
        }

        // fixed nodes
        for( auto& node_i : mListOfFixedNodes )
        {
            // Get node information
            int i = node_i->GetValue( MAPPING_ID );

            int counter_col = 0;
            for( auto& node_j : mListOfActiveNodes )
            {
                int j = node_j->GetValue( MAPPING_ID );
                double val = mMappingMatrix(i,j);
                if (val != 0.0)
                    modifiedMappingMatrix.insert_element( counter_row, counter_col, penalty_factor*val );
                counter_col ++;
            }

            x_variables_modified[counter_row] = 0;
            y_variables_modified[counter_row] = 0;
            z_variables_modified[counter_row] = 0;

            if (!mFix_in_X)
                x_variables_modified[counter_row] = penalty_factor * x_variables_in_geometry_space[i];
            if (!mFix_in_Y)
                y_variables_modified[counter_row] = penalty_factor * y_variables_in_geometry_space[i];
            if (!mFix_in_Z)
                z_variables_modified[counter_row] = penalty_factor * z_variables_in_geometry_space[i];

            counter_row ++;
        }

        std::cout << "> Time needed for assembling modifiedMappingMatrix: " << assembling_time.elapsed() << " s" << std::endl;

        // write mapping matrix

        // writeMatrixToFile( modifiedMappingMatrix,
        //                    numberOfRowsInModifiedSystem,
        //                    numberOfColsInModifiedSystem,
        //                    "modifiedMappingMatrixReduced.txt");

        // writeVectorToFile ( x_variables_modified,
        //                     "x_variables_modified_red.txt");
        // writeVectorToFile ( y_variables_modified,
        //                     "y_variables_modified_red.txt");
        // writeVectorToFile ( z_variables_modified,
        //                     "z_variables_modified_red.txt");

        // writeMatrixToFile( mMappingMatrix,
        //                    mNumberOfDesignVariables,
        //                    mNumberOfDesignVariables,
        //                    "originalMappingMatrix.txt");

        boost::timer solving_time;
        std::cout << "\n> Starting to solve modified System..." << std::endl;
        Vector x_variables_in_design_space_affected, y_variables_in_design_space_affected, z_variables_in_design_space_affected;

        x_variables_in_design_space_affected.resize(mListOfActiveNodes.size(),0.0);
        y_variables_in_design_space_affected.resize(mListOfActiveNodes.size(),0.0);
        z_variables_in_design_space_affected.resize(mListOfActiveNodes.size(),0.0);

        linear_solver->Solve(modifiedMappingMatrix, x_variables_in_design_space_affected, x_variables_modified);
        linear_solver->Solve(modifiedMappingMatrix, y_variables_in_design_space_affected, y_variables_modified);
        linear_solver->Solve(modifiedMappingMatrix, z_variables_in_design_space_affected, z_variables_modified);

        std::cout << "> Time needed for solving modified System: " << solving_time.elapsed() << " s" << std::endl;

        boost::timer backassambling_time;
        std::cout << "\n> Starting backward-assembling..." << std::endl;

        // backward-assembling in the design-vector
        counter_row = 0;
        for( auto& node_i : mListOfActiveNodes )
        {
            int i = node_i->GetValue( MAPPING_ID );
            x_variables_in_design_space[i] = x_variables_in_design_space_affected[counter_row];
            y_variables_in_design_space[i] = y_variables_in_design_space_affected[counter_row];
            z_variables_in_design_space[i] = z_variables_in_design_space_affected[counter_row];

            counter_row ++;
        }

        // writeVectorToFile ( x_variables_in_design_space,
        //                     "x_variables_in_design_space_red.txt");
        // writeVectorToFile ( y_variables_in_design_space,
        //                     "y_variables_in_design_space_red.txt");
        // writeVectorToFile ( z_variables_in_design_space,
        //                     "z_variables_in_design_space_red.txt");

        std::cout << "> Time needed for backward-assembling: " << backassambling_time.elapsed() << " s" << std::endl;

    }

    // --------------------------------------------------------------------------
    void writeMatrixToFile( Matrix matrix, int nrow, int ncol, const char* filename )
    {
        std::ofstream file (filename);
        if (file.is_open())
        {
            for (int i=0;i<nrow;i++)
            {
                for (int j=0;j<ncol;j++)
                {
                    double val = matrix(i,j);
                    file << val << ", ";
                }
                file << "\n";
            }
            file.close();
        }
    }

    // --------------------------------------------------------------------------
    void writeVectorToFile( NodeVector vec, const char* filename )
    {
        std::ofstream file (filename);
        if (file.is_open())
        {
            for (auto& node : vec )
            {
                int val = node->GetValue(MAPPING_ID);
                file << val << " \n";
            }
            file.close();
        }
    }

    // --------------------------------------------------------------------------
    void writeVectorToFile( Vector vec, const char* filename )
    {

        std::ofstream file (filename);
        if (file.is_open())
        {
            for ( unsigned int i=0;i<vec.size();i++)
            {
                double val = vec[i];
                file << val << " \n";
            }
            file.close();
        }
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
        return "MapperVertexMorphingRigidBody";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperVertexMorphingRigidBody";
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
    FilterFunction::Pointer mpFilterFunction;

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
    std::string mFilterType;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;
    bool mConsistentBackwardMapping;
    NodeVector mListOfRigidNodes;
    NodeVector mListOfFixedNodes;
    unsigned int mNumberOfRigidNodes = 0;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesOfDesignSurface;
    KDTree::Pointer mpSearchTree;

    // ==============================================================================
    // Variables for mapping
    // ==============================================================================
    CompressedMatrixType mMappingMatrix;
    Vector x_variables_in_design_space, y_variables_in_design_space, z_variables_in_design_space;
    Vector x_variables_in_geometry_space, y_variables_in_geometry_space, z_variables_in_geometry_space;
    Vector mTranslationVector;
    Matrix mRotationMatrix;
    NodeVector mListOfActiveNodes;
    double mControlSum = 0.0;

    boost::python::list mrRigidRegions;
    boost::python::list mrFixedRegions;
    bool mFix_in_X, mFix_in_Y, mFix_in_Z;
    Parameters mMapperSettings;


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
//      MapperVertexMorphingRigidBody& operator=(MapperVertexMorphingRigidBody const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingRigidBody(MapperVertexMorphingRigidBody const& rOther);


    ///@}

}; // Class MapperVertexMorphingRigidBody

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_H
