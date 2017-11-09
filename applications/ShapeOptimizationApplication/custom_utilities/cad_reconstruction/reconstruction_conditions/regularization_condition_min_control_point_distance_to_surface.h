// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef REGULARIZATION_CONDITION_MIN_CP_DISTANCE_TO_SURFACE_H
#define REGULARIZATION_CONDITION_MIN_CP_DISTANCE_TO_SURFACE_H

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
#include "regularization_condition_base.h"

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

class MinimalControlPointDistanceToSurfaceCondition : public RegularizationCondition
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Pointer definition of MinimalControlPointDistanceToSurfaceCondition
    KRATOS_CLASS_POINTER_DEFINITION(MinimalControlPointDistanceToSurfaceCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MinimalControlPointDistanceToSurfaceCondition( ReconstructionDataBase& data_base, Parameters& reconstruction_parameters )
    : mrReconstructionDataBase( data_base ),
      mAlphaValue( reconstruction_parameters["solution_parameters"]["regularization_parameters"]["alpha"].GetDouble() )
    {
    }

    /// Destructor.
    virtual ~MinimalControlPointDistanceToSurfaceCondition()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================  
    void Initialize()
    {
    }  
    
    // --------------------------------------------------------------------------
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS )
    {
        std::cout << "> Starting to compute LHS contribution of alpha-regularization (minimal control point distance to surface)..." << std::endl;  
        
        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            std::vector<double> greville_abscissae_in_u_direction;
            std::vector<double> greville_abscissae_in_v_direction;
    
            patch_i.ComputeGrevilleAbscissae( greville_abscissae_in_u_direction, greville_abscissae_in_v_direction );
    
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
            {
                int control_point_vector_index = &control_point_i-&(patch_i.GetSurfaceControlPoints()[0]);
                if(control_point_i.IsRelevantForReconstruction())
                {
                    int control_point_id = control_point_i.GetEquationId();

                    array_1d<double,2> greville_parameters;
                    greville_parameters[0] = greville_abscissae_in_u_direction[control_point_vector_index];
                    greville_parameters[1] = greville_abscissae_in_v_direction[control_point_vector_index];
                    array_1d<int,2> parameter_spans = patch_i.ComputeSurfaceKnotSpans( greville_parameters );
                    
                    std::vector<double> nurbs_function_values = patch_i.EvaluateNURBSFunctions( parameter_spans, greville_parameters );
                    std::vector<ControlPoint*> cps_affecting_point_on_greville_abscissae = patch_i.GetPointersToAffectedControlPoints( parameter_spans, greville_parameters );
                    int number_of_control_points_affecting_greville_point = cps_affecting_point_on_greville_abscissae.size();

                    // LHS contribution 1 (influence of control point under consideration )
                    LHS(3*control_point_id+0, 3*control_point_id+0) += mAlphaValue;
                    LHS(3*control_point_id+1, 3*control_point_id+1) += mAlphaValue;
                    LHS(3*control_point_id+2, 3*control_point_id+2) += mAlphaValue;

                    for(int row_itr=0; row_itr<number_of_control_points_affecting_greville_point; row_itr++)
                    {
                        ControlPoint* row_control_point = cps_affecting_point_on_greville_abscissae[row_itr];
                        if(row_control_point->IsRelevantForReconstruction())
                        {                        
                            int row_id = row_control_point->GetEquationId();
                            double R_row = nurbs_function_values[row_itr];

                            // LHS contribution 2 (mixed influence of control point under consideration & control points affecting point on greville abscissae)
                            LHS(3*control_point_id+0, 3*row_id+0) -= mAlphaValue * R_row;
                            LHS(3*control_point_id+1, 3*row_id+1) -= mAlphaValue * R_row;
                            LHS(3*control_point_id+2, 3*row_id+2) -= mAlphaValue * R_row;
                            // LHS contribution 3 (mixed influence of control point under consideration & control points affecting point on greville abscissae)
                            LHS(3*row_id+0, 3*control_point_id+0) -= mAlphaValue * R_row;
                            LHS(3*row_id+1, 3*control_point_id+1) -= mAlphaValue * R_row;
                            LHS(3*row_id+2, 3*control_point_id+2) -= mAlphaValue * R_row;                            
                
                            for(int collumn_itr=0; collumn_itr<number_of_control_points_affecting_greville_point; collumn_itr++)
                            {      
                                ControlPoint* collumn_control_point = cps_affecting_point_on_greville_abscissae[collumn_itr];
                                if(collumn_control_point->IsRelevantForReconstruction())
                                {    
                                    int collumn_id = collumn_control_point->GetEquationId();
                                    double R_collumn = nurbs_function_values[collumn_itr];

                                    // LHS contribution 4 (influence of control points affecting point on greville abscissae)
                                    LHS(3*row_id+0, 3*collumn_id+0) += mAlphaValue * R_row*R_collumn;
                                    LHS(3*row_id+1, 3*collumn_id+1) += mAlphaValue * R_row*R_collumn;
                                    LHS(3*row_id+2, 3*collumn_id+2) += mAlphaValue * R_row*R_collumn;
                                }
                            }
                        }                   
                    }
                }
            }
        }
        std::cout << "> Finished computing LHS contribution of alpha-regularization (minimal control point distance to surface)..." << std::endl;          
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS )
    {        
        std::cout << "> Starting to compute RHS contribution of alpha-regularization (minimal control point distance to surface)..." << std::endl;  
    
        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            std::vector<double> greville_abscissae_in_u_direction;
            std::vector<double> greville_abscissae_in_v_direction;
    
            patch_i.ComputeGrevilleAbscissae( greville_abscissae_in_u_direction, greville_abscissae_in_v_direction );
    
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
            {
                int control_point_vector_index = &control_point_i-&(patch_i.GetSurfaceControlPoints()[0]);
                if(control_point_i.IsRelevantForReconstruction())
                {
                    int control_point_id = control_point_i.GetEquationId();

                    array_1d<double,2> greville_parameters;
                    greville_parameters[0] = greville_abscissae_in_u_direction[control_point_vector_index];
                    greville_parameters[1] = greville_abscissae_in_v_direction[control_point_vector_index];
                    array_1d<int,2> parameter_spans = patch_i.ComputeSurfaceKnotSpans( greville_parameters );

                    Point<3> greville_point;
                    patch_i.EvaluateSurfacePoint( greville_parameters, greville_point );                    
                    
                    std::vector<double> nurbs_function_values = patch_i.EvaluateNURBSFunctions( parameter_spans, greville_parameters );
                    std::vector<ControlPoint*> cps_affecting_point_on_greville_abscissae = patch_i.GetPointersToAffectedControlPoints( parameter_spans, greville_parameters );
                    int number_of_control_points_affecting_greville_point = cps_affecting_point_on_greville_abscissae.size();

                    // RHS contribution 1 (influence of control point under consideration )
                    RHS(3*control_point_id+0) -= mAlphaValue * (control_point_i.GetX() - greville_point.X());
                    RHS(3*control_point_id+1) -= mAlphaValue * (control_point_i.GetY() - greville_point.Y());
                    RHS(3*control_point_id+2) -= mAlphaValue * (control_point_i.GetZ() - greville_point.Z());

                    for(int row_itr=0; row_itr<number_of_control_points_affecting_greville_point; row_itr++)
                    {
                        ControlPoint* row_control_point = cps_affecting_point_on_greville_abscissae[row_itr];
                        if(row_control_point->IsRelevantForReconstruction())
                        {                        
                            int row_id = row_control_point->GetEquationId();
                            double R_row = nurbs_function_values[row_itr];

                            // RHS contribution 2 (influence of control points affecting point on greville abscissae)
                            RHS(3*row_id+0) += mAlphaValue * R_row * (control_point_i.GetX() - greville_point.X());
                            RHS(3*row_id+1) += mAlphaValue * R_row * (control_point_i.GetY() - greville_point.Y());
                            RHS(3*row_id+2) += mAlphaValue * R_row * (control_point_i.GetZ() - greville_point.Z());                          
                        }                   
                    }
                }
            }
        } 
        std::cout << "> Finished computing RHS contribution of alpha-regularization (minimal control point distance to surface)..." << std::endl;                  
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
        return "MinimalControlPointDistanceToSurfaceCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "MinimalControlPointDistanceToSurfaceCondition";
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ReconstructionDataBase& mrReconstructionDataBase;
    double mAlphaValue;
    std::string mSolutionStrategy;

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
    //      MinimalControlPointDistanceToSurfaceCondition& operator=(MinimalControlPointDistanceToSurfaceCondition const& rOther);

    /// Copy constructor.
    //      MinimalControlPointDistanceToSurfaceCondition(MinimalControlPointDistanceToSurfaceCondition const& rOther);

    ///@}

}; // Class MinimalControlPointDistanceToSurfaceCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // REGULARIZATION_CONDITION_MIN_CP_DISTANCE_TO_SURFACE_H
