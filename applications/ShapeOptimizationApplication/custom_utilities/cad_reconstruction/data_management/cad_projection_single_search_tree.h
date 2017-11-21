// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_PROJECTION_SINGLE_SEARCH_TREE_H
#define CAD_PROJECTION_SINGLE_SEARCH_TREE_H

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
class CADProjectionSingleSearchTree : public CADProjectionBase
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

    /// Pointer definition of CADProjectionSingleSearchTree
    KRATOS_CLASS_POINTER_DEFINITION(CADProjectionSingleSearchTree);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADProjectionSingleSearchTree( PatchVector& patch_vector, Parameters projection_parameters )
    : mrPatchVector( patch_vector ),
      mProjectionParameters( projection_parameters )
    {      
    }

    /// Destructor.
    virtual ~CADProjectionSingleSearchTree()
    {
    }

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        std::cout << "\n> Initializing CAD projection..." << std::endl;           
        boost::timer timer;

        bool initialize_using_greville_abscissae = mProjectionParameters["automatic_initialization_using_greville_abscissae"].GetBool();
        if(initialize_using_greville_abscissae)
          CreateCADPointCloudBasedOnGrevilleAbscissae();
        else
          CreateCADPointCloudBasedOnManualInput();
        CreateSearchTreeWithCADPointCloud();

    std::cout << "> Time needed initializing CAD projection: " << timer.elapsed() << " s" << std::endl;    
    }  

    // --------------------------------------------------------------------------
    void DetermineNearestCADPoint( NodeType& PointOfInterest,
                                   array_1d<double,2>& parameter_values_of_nearest_point,
                                   int& patch_index_of_nearest_point ) override
    {
        // 1) Coarse search in the point cloud
        NodeType::Pointer nearest_point = mpSearchTree->SearchNearestPoint( PointOfInterest );

        // 2) Detailed projection using Newton-Raphson
        
        // Recover CAD information from lists representing point cloud
        parameter_values_of_nearest_point = mParameterValuesOfNodesInGlobalSearchTree[nearest_point->Id()-1];
        patch_index_of_nearest_point = mPatchIndicesOfNodesInGlobalSearchTree[nearest_point->Id()-1];
        Patch& patch_of_nearest_point = mrPatchVector[patch_index_of_nearest_point];
        
        OptimizeGuessWithNewtonRaphson( PointOfInterest, 
                                        *nearest_point, 
                                        parameter_values_of_nearest_point, 
                                        patch_of_nearest_point,
                                        mProjectionParameters );
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
    return "CADProjectionSingleSearchTree";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
    rOStream << "CADProjectionSingleSearchTree";
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

    // Variables for spatial search
    unsigned int mNumberOfNodesInGlobalSearchTree = 0;
    NodeVector mNodesForGlobalSearchTree;
    std::vector<array_1d<double,2>> mParameterValuesOfNodesInGlobalSearchTree; 
    std::vector<int> mPatchIndicesOfNodesInGlobalSearchTree;
    unsigned int mBucketSize = 100;
    KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private member functions
    ///@{   
    void CreateCADPointCloudBasedOnGrevilleAbscissae()
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
            
            // Points of (refined) Greville abscissae are rendered into a point cloud
            number_of_greville_points = greville_abscissae_in_u_direction.size();
            array_1d<double,2> point_in_parameter_space;
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
                    mParameterValuesOfNodesInGlobalSearchTree.push_back(point_in_parameter_space);
                    mPatchIndicesOfNodesInGlobalSearchTree.push_back(index_in_patch_vector);
                }
            }
        }      
    }    
    
    // --------------------------------------------------------------------------
    void CreateCADPointCloudBasedOnManualInput()
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
                        mParameterValuesOfNodesInGlobalSearchTree.push_back(point_in_parameter_space);
                        mPatchIndicesOfNodesInGlobalSearchTree.push_back(index_in_patch_vector);
                    }
                }
            }
        }      
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithCADPointCloud()
    {
      std::cout << "> Creating search tree..." << std::endl;           
      mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mNodesForGlobalSearchTree.begin(), mNodesForGlobalSearchTree.end(), mBucketSize));
      std::cout << "> Search tree created using " << mNodesForGlobalSearchTree.size() << " points on the CAD model." << std::endl;                 
    }       

    /// Assignment operator.
    //      CADProjectionSingleSearchTree& operator=(CADProjectionSingleSearchTree const& rOther);

    /// Copy constructor.
    //      CADProjectionSingleSearchTree(CADProjectionSingleSearchTree const& rOther);

}; // Class CADProjectionSingleSearchTree
} // namespace Kratos.

#endif // CAD_PROJECTION_SINGLE_SEARCH_TREE_H
