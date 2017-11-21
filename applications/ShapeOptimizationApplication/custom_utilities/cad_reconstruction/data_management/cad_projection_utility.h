// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_PROJECTION_UTILITY_H
#define CAD_PROJECTION_UTILITY_H

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
#include "../basic_nurbs_brep_handling/patch.h"
#include "cad_projection_base.h"
#include "cad_projection_single_search_tree.h"
#include "cad_projection_multiple_search_trees.h"

// ==============================================================================

namespace Kratos
{
class CADProjectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef std::vector<Patch> PatchVector; 

    /// Pointer definition of CADProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(CADProjectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADProjectionUtility( PatchVector& patch_vector, Parameters projection_parameters )
    : mrPatchVector( patch_vector )
    {
        if(projection_parameters["projection_strategy"].GetString().compare("single_search_tree")==0)
            mpProjectionStrategy = CADProjectionBase::Pointer (new CADProjectionSingleSearchTree( patch_vector, projection_parameters ));
        else if(projection_parameters["projection_strategy"].GetString().compare("multiple_search_tree")==0)
            mpProjectionStrategy = CADProjectionBase::Pointer (new CADProjectionMultipleSearchTrees( patch_vector, projection_parameters ));
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Specified projection strategy not implemented!", "");
    }

    /// Destructor.
    virtual ~CADProjectionUtility()
    {
    }

    // --------------------------------------------------------------------------
    void Initialize()
    {
        mpProjectionStrategy->Initialize();
    }

    // --------------------------------------------------------------------------
    void DetermineNearestCADPoint( NodeType& PointOfInterest, array_1d<double,2>& parameter_values_of_nearest_point, int& patch_index_of_nearest_point )
    {
        mpProjectionStrategy->DetermineNearestCADPoint(PointOfInterest, parameter_values_of_nearest_point, patch_index_of_nearest_point );
    }
    
    // --------------------------------------------------------------------------
    void DetermineNearestCADPointInGeometrySpace( NodeType& PointOfInterest, array_1d<double,3>& coordinates_of_nearest_cad_point )
    {
        array_1d<double,2> parameter_values_of_nearest_point;
        int patch_index_of_nearest_point = -1;

        DetermineNearestCADPoint( PointOfInterest, parameter_values_of_nearest_point, patch_index_of_nearest_point );

        Point<3> nearest_cad_point;
        mrPatchVector[patch_index_of_nearest_point].EvaluateSurfacePoint( parameter_values_of_nearest_point, nearest_cad_point );

        coordinates_of_nearest_cad_point[0] = nearest_cad_point[0];
        coordinates_of_nearest_cad_point[1] = nearest_cad_point[1];
        coordinates_of_nearest_cad_point[2] = nearest_cad_point[2];      
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
    return "CADProjectionUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
    rOStream << "CADProjectionUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }
    
    // ==============================================================================
  private:

    // Variables initialized by constructor
    CADProjectionBase::Pointer mpProjectionStrategy;
    PatchVector& mrPatchVector;

    /// Assignment operator.
    //      CADProjectionUtility& operator=(CADProjectionUtility const& rOther);

    /// Copy constructor.
    //      CADProjectionUtility(CADProjectionUtility const& rOther);

}; // Class CADProjectionUtility
} // namespace Kratos.

#endif // CAD_PROJECTION_UTILITY_H
