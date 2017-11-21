// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_PROJECTION_MULTIPLE_SEARCH_TREES_H
#define CAD_PROJECTION_MULTIPLE_SEARCH_TREES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "reconstruction_data_base.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"
#include "cad_projection_base.h"

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
class CADProjectionMultipleSearchTrees : public CADProjectionBase
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;    
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;    
    typedef std::vector<Patch> PatchVector; 
    typedef std::vector<double> DoubleVector;
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;    

    /// Pointer definition of CADProjectionMultipleSearchTrees
    KRATOS_CLASS_POINTER_DEFINITION(CADProjectionMultipleSearchTrees);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADProjectionMultipleSearchTrees( PatchVector& patch_vector, Parameters projection_parameters )
    : mrPatchVector( patch_vector ),
      mProjectionParameters( projection_parameters )
    {
        mSearchRadius = mProjectionParameters["search_radius_with_multiple_trees"].GetDouble();
    }

    /// Destructor.
    virtual ~CADProjectionMultipleSearchTrees()
    {
    }

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        std::cout << "\n> Initializing CAD projection..." << std::endl;           
        boost::timer timer;

        bool initialize_using_greville_abscissae = mProjectionParameters["automatic_initialization_using_greville_abscissae"].GetBool();
        if(initialize_using_greville_abscissae)
          CreateCADPointCloudsBasedOnGrevilleAbscissae();
        else
          CreateCADPointsCloudBasedOnManualInput();
        CreateSearchTreesForEachPointCloud();

    std::cout << "> Time needed initializing CAD projection: " << timer.elapsed() << " s" << std::endl;    
    }  

    // --------------------------------------------------------------------------
    void DetermineNearestCADPoint( NodeType& PointOfInterest,
                                   array_1d<double,2>& parameter_values_of_nearest_point,
                                   int& patch_index_of_nearest_point ) override
    {
        NodeType::Pointer nearest_point_in_global_tree = mpGlobalSearchTree->SearchNearestPoint( PointOfInterest );

        parameter_values_of_nearest_point = mParameterValuesOfNodesInGlobalSearchTree[nearest_point_in_global_tree->Id()-1];
        patch_index_of_nearest_point = mPatchIndicesOfNodesInGlobalSearchTree[nearest_point_in_global_tree->Id()-1];
        Patch& patch_of_nearest_point_in_global_tree = mrPatchVector[patch_index_of_nearest_point];

        OptimizeGuessWithNewtonRaphson( PointOfInterest, 
                                        *nearest_point_in_global_tree, 
                                        parameter_values_of_nearest_point, 
                                        patch_of_nearest_point_in_global_tree,
                                        mProjectionParameters );
       
        // If nearest point is inside --> no search check necessary
        if(patch_of_nearest_point_in_global_tree.IsPointInside(parameter_values_of_nearest_point))
            return;
        // If nearest point is outside, we loop over all patch trees and look for the nearest point until one is inside
        // Note that only nearest points inside a certain search readius will be further analyzed by NewtonRaphson (this is to prevent unecessary analysis for patches far away) 
        else
        {
          for(int tree_index = 0; tree_index < mListOfLocalSearchTrees.size(); tree_index++)
          {
              NodeType::Pointer nearest_point_in_local_tree = mListOfLocalSearchTrees[tree_index]->SearchNearestPoint( PointOfInterest );
              parameter_values_of_nearest_point = mParameterValuesOfNodesInGlobalSearchTree[nearest_point_in_local_tree->Id()-1];
              patch_index_of_nearest_point = mPatchIndicesOfNodesInGlobalSearchTree[nearest_point_in_local_tree->Id()-1];
              Patch& patch_of_nearest_point_in_local_tree = mrPatchVector[patch_index_of_nearest_point];

            if(DistanceBetweenNodes(PointOfInterest, *nearest_point_in_local_tree) < mSearchRadius)
            {
              OptimizeGuessWithNewtonRaphson( PointOfInterest, 
                                              *nearest_point_in_local_tree, 
                                              parameter_values_of_nearest_point, 
                                              patch_of_nearest_point_in_local_tree, 
                                              mProjectionParameters );
              if(patch_of_nearest_point_in_local_tree.IsPointInside(parameter_values_of_nearest_point))
                return;
            }
          }
        }
        KRATOS_THROW_ERROR(std::runtime_error, "Unable to find CAD nearest neighbour inside trimming boundary", "")
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
    return "CADProjectionMultipleSearchTrees";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
    rOStream << "CADProjectionMultipleSearchTrees";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

    // ==============================================================================
  private:

    
    // Variables initialized by constructor
    PatchVector& mrPatchVector;
    Parameters mProjectionParameters;
    double mSearchRadius;

    // Variables for spatial search
    unsigned int mNumberOfNodesInGlobalSearchTree = 0;
    NodeVector mNodesForGlobalSearchTree;
    std::vector< NodeVector > mNodesForLocalPatchTree;
    std::vector<array_1d<double,2>> mParameterValuesOfNodesInGlobalSearchTree; 
    std::vector<int> mPatchIndicesOfNodesInGlobalSearchTree;
    unsigned int mBucketSize = 100;
    KDTree::Pointer mpGlobalSearchTree;
    std::vector< KDTree::Pointer > mListOfLocalSearchTrees;

    ///@}
    ///@name Private member functions
    ///@{   
    void CreateCADPointCloudsBasedOnGrevilleAbscissae()
    {
      for (auto & patch_i : mrPatchVector)
      {
            int index_in_patch_vector = &patch_i - &mrPatchVector[0];

            // Computation of greville points
            std::vector<double> greville_abscissae_in_u_direction;
            std::vector<double> greville_abscissae_in_v_direction;
            patch_i.ComputeGrevilleAbscissae( greville_abscissae_in_u_direction, greville_abscissae_in_v_direction );

            // Greville abscissae is refined if specified
            int number_of_refinements = mProjectionParameters["refinement_iterations_of_greville_abscissae"].GetInt();      
            int number_of_greville_points = greville_abscissae_in_u_direction.size();
            std::vector<double> refined_greville_abscissae_in_u_direction;
            std::vector<double> refined_greville_abscissae_in_v_direction;            
            for(int refinement_itr=0; refinement_itr<number_of_refinements; refinement_itr++)
            {
                patch_i.RefineGrevilleAbscissae( greville_abscissae_in_u_direction, 
                                                 greville_abscissae_in_v_direction,
                                                 refined_greville_abscissae_in_u_direction, 
                                                 refined_greville_abscissae_in_v_direction );
                
                greville_abscissae_in_u_direction = refined_greville_abscissae_in_u_direction;
                greville_abscissae_in_v_direction = refined_greville_abscissae_in_v_direction;
            }
            
            // Points of (refined) Greville abscissae are rendered into global and local point cloud
            number_of_greville_points = greville_abscissae_in_u_direction.size();
            array_1d<double,2> point_in_parameter_space;

            NodeVector empty_list_of_nodes;
            mNodesForLocalPatchTree.push_back(empty_list_of_nodes);  

            for(int i=0; i<number_of_greville_points; i++)
            {
                point_in_parameter_space[0] = greville_abscissae_in_u_direction[i];
                point_in_parameter_space[1] = greville_abscissae_in_v_direction[i];

                bool point_is_inside = patch_i.IsPointInside(point_in_parameter_space);
                if(point_is_inside)
                {
                    ++mNumberOfNodesInGlobalSearchTree;					
                    Point<3> cad_point_coordinates;
                    patch_i.EvaluateSurfacePoint( point_in_parameter_space, cad_point_coordinates );

                    NodeType::Pointer new_cad_node = Node <3>::Pointer(new Node<3>(mNumberOfNodesInGlobalSearchTree, cad_point_coordinates));

                    mNodesForGlobalSearchTree.push_back(new_cad_node);
                    mNodesForLocalPatchTree[index_in_patch_vector].push_back(new_cad_node);
                    mParameterValuesOfNodesInGlobalSearchTree.push_back(point_in_parameter_space);
                    mPatchIndicesOfNodesInGlobalSearchTree.push_back(index_in_patch_vector);
                }
            }
        }      
    }    
    
    // --------------------------------------------------------------------------
    void CreateCADPointsCloudBasedOnManualInput()
    {
      int u_resolution = mProjectionParameters["parameter_resolution_for_manual_initialization"][0].GetInt();
      int v_resolution =  mProjectionParameters["parameter_resolution_for_manual_initialization"][1].GetInt();

      for (auto & patch_i : mrPatchVector)
      {
            int index_in_patch_vector = &patch_i - &mrPatchVector[0];
            DoubleVector& knot_vec_u_i = patch_i.GetSurfaceKnotVectorU();
            DoubleVector& knot_vec_v_i = patch_i.GetSurfaceKnotVectorV();
            std::cout << "> Processing Patch with brep_id " << patch_i.GetId() << std::endl;
      
            double u_min = knot_vec_u_i[0];
            double u_max = knot_vec_u_i[knot_vec_u_i.size()-1];
            double v_min = knot_vec_v_i[0];
            double v_max = knot_vec_v_i[knot_vec_v_i.size()-1];
            double delta_u = (u_max-u_min) / u_resolution;
            double delta_v = (v_max-v_min) / v_resolution;

            // Loop over all u & v according to specified resolution
            array_1d<double,2> point_in_parameter_space;

            NodeVector empty_list_of_nodes;
            mNodesForLocalPatchTree.push_back(empty_list_of_nodes);  

            for(int i=0; i<=u_resolution; i++)
            {
                point_in_parameter_space[0] = u_min + i*delta_u;

                for(int j=0; j<=v_resolution; j++)
                {
                    point_in_parameter_space[1] = v_min + j*delta_v;

                    bool point_is_inside = patch_i.IsPointInside(point_in_parameter_space);
                    if(point_is_inside)
                    {
                        ++mNumberOfNodesInGlobalSearchTree;					
                        Point<3> cad_point_coordinates;
                        patch_i.EvaluateSurfacePoint( point_in_parameter_space, cad_point_coordinates );

                         NodeType::Pointer new_cad_node = Node <3>::Pointer(new Node<3>(mNumberOfNodesInGlobalSearchTree, cad_point_coordinates));

                        mNodesForGlobalSearchTree.push_back(new_cad_node);
                        mNodesForLocalPatchTree[index_in_patch_vector].push_back(new_cad_node);                        
                        mParameterValuesOfNodesInGlobalSearchTree.push_back(point_in_parameter_space);
                        mPatchIndicesOfNodesInGlobalSearchTree.push_back(index_in_patch_vector);
                    }
                }
            }
        }      
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreesForEachPointCloud()
    {
        std::cout << "> Creating multiple search trees..." << std::endl;
        
        // Creating global search tree
        mpGlobalSearchTree = boost::shared_ptr<KDTree>(new KDTree(mNodesForGlobalSearchTree.begin(), mNodesForGlobalSearchTree.end(), mBucketSize));   

        // Creating local search tree
        for(auto & list_of_nodes_in_local_patch : mNodesForLocalPatchTree)
            mListOfLocalSearchTrees.push_back( boost::shared_ptr<KDTree>(new KDTree(list_of_nodes_in_local_patch.begin(), list_of_nodes_in_local_patch.end(), mBucketSize)) );

        std::cout << "> Search trees created." << std::endl;                 
    }
    
    // --------------------------------------------------------------------------
    double DistanceBetweenNodes(NodeType& node_1, NodeType& node_2)
    {
      return std::sqrt( (node_1.X() - node_2.X()) * (node_1.X() - node_2.X()) +
                        (node_1.Y() - node_2.Y()) * (node_1.Y() - node_2.Y()) + 
                        (node_1.Z() - node_2.Z()) * (node_1.Z() - node_2.Z()) );
    }         

    /// Assignment operator.
    //      CADProjectionMultipleSearchTrees& operator=(CADProjectionMultipleSearchTrees const& rOther);

    /// Copy constructor.
    //      CADProjectionMultipleSearchTrees(CADProjectionMultipleSearchTrees const& rOther);

}; // Class CADProjectionMultipleSearchTrees
} // namespace Kratos.

#endif // CAD_PROJECTION_MULTIPLE_SEARCH_TREE_H
