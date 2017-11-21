// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONDITION_DISTANCE_MINIMIZATION_H
#define RECONSTRUCTION_CONDITION_DISTANCE_MINIMIZATION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "reconstruction_condition_base.h"

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

class DistanceMinimizationCondition : public ReconstructionCondition
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;        

    /// Pointer definition of DistanceMinimizationCondition
    KRATOS_CLASS_POINTER_DEFINITION(DistanceMinimizationCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceMinimizationCondition( NodeType&  fem_node,
                                   Patch& patch,
                                   array_1d<double,2> param_values )
    : mrFemNode( fem_node ),
      mrAffectedPatch( patch ),
      mParmeterValues( param_values )
    {
        mParameterSpans = mrAffectedPatch.ComputeSurfaceKnotSpans( mParmeterValues );
    }

    /// Destructor.
    virtual ~DistanceMinimizationCondition()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void FlagControlPointsRelevantForReconstruction()
    {
        mrAffectedPatch.FlagAffectedControlPointsForReconstruction(  mParameterSpans, mParmeterValues );
    }    
    
    // --------------------------------------------------------------------------
    void Initialize()
    {        
        mNurbsFunctionValues = mrAffectedPatch.EvaluateNURBSFunctions( mParameterSpans, mParmeterValues );
        mAffectedControlPoints = mrAffectedPatch.GetPointersToAffectedControlPoints( mParameterSpans, mParmeterValues );
        mEquationIdsOfAffectedControlPoints = mrAffectedPatch.GetEquationIdsOfAffectedControlPoints( mParameterSpans, mParmeterValues );
        mNumberOfLocalEquationIds = mEquationIdsOfAffectedControlPoints.size();         
    }    
    
    // --------------------------------------------------------------------------
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) override
    {           
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIds; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPoints[row_itr];
            double R_row = mNurbsFunctionValues[row_itr];

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIds; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPoints[collumn_itr];
                double R_collumn = mNurbsFunctionValues[collumn_itr];

                LHS( 3*row_id+0, 3*collumn_id+0 ) += R_row * R_collumn;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += R_row * R_collumn;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += R_row * R_collumn;
            }
        }
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        // Compute RHS
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIds; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPoints[row_itr];
            double R_row = mNurbsFunctionValues[row_itr];

            // Computation of -R*\hat{P_F}
            Vector fem_contribution = ZeroVector(3);
            array_1d<double,3>& node_disp = mrFemNode.FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE);
            fem_contribution(0) = R_row * ( mrFemNode.X() + node_disp[0] );
            fem_contribution(1) = R_row * ( mrFemNode.Y() + node_disp[1] );
            fem_contribution(2) = R_row * ( mrFemNode.Z() + node_disp[2] );
            
            // Computation of RR*\hat{P_C}
            Vector cad_cad_contribution = ZeroVector(3);                                            
            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIds; collumn_itr++)
            {                
                double R_collumn = mNurbsFunctionValues[collumn_itr];

                cad_cad_contribution(0) += R_row * R_collumn * mAffectedControlPoints[collumn_itr]->GetX();
                cad_cad_contribution(1) += R_row * R_collumn * mAffectedControlPoints[collumn_itr]->GetY();
                cad_cad_contribution(2) += R_row * R_collumn * mAffectedControlPoints[collumn_itr]->GetZ();
            }  

            // Computation of complete RHS contribution: rhs_contribution = -( R*\hat{P_C} - \hat{P_F})R
            RHS( 3*row_id+0 ) -= ( cad_cad_contribution[0] - fem_contribution[0]);
            RHS( 3*row_id+1 ) -= ( cad_cad_contribution[1] - fem_contribution[1]);
            RHS( 3*row_id+2 ) -= ( cad_cad_contribution[2] - fem_contribution[2]);             
        }
    }

    // --------------------------------------------------------------------------
    void DetermineFECoordinatesInUndeformedConfiguration( array_1d<double,3>& fe_point_coords ) override
    {
        fe_point_coords[0] = mrFemNode.X();
        fe_point_coords[1] = mrFemNode.Y();  
        fe_point_coords[2] = mrFemNode.Z();                        
    }

    // --------------------------------------------------------------------------
    void DetermineFECoordinatesInDeformedConfiguration( Variable<array_1d<double,3>> shape_change_variable, array_1d<double,3>& fe_point_coords_in_deformed_configuration ) override
    {
        array_1d<double,3>& displacement = mrFemNode.FastGetSolutionStepValue(shape_change_variable);
        fe_point_coords_in_deformed_configuration[0] = mrFemNode.X() + displacement[0];
        fe_point_coords_in_deformed_configuration[1] = mrFemNode.Y() + displacement[1];
        fe_point_coords_in_deformed_configuration[2] = mrFemNode.Z() + displacement[2];        
    }
    // --------------------------------------------------------------------------
    void DetermineCADCoordinatesInUndeformedConfiguration(array_1d<double,3>& cad_point_coords) override
    {
        Point<3> cad_point_in_deformed_configuration;
        mrAffectedPatch.EvaluateSurfacePoint(mParmeterValues, cad_point_in_deformed_configuration);

        array_1d<double,3> displacement;
        mrAffectedPatch.EvaluateSurfaceDisplacement(mParmeterValues, displacement);

        cad_point_coords[0] = cad_point_in_deformed_configuration(0) - displacement[0];
        cad_point_coords[1] = cad_point_in_deformed_configuration(1) - displacement[1];
        cad_point_coords[2] = cad_point_in_deformed_configuration(2) - displacement[2];        
    }
    // --------------------------------------------------------------------------
    void DetermineCADCoordinatesInDeformedConfiguration(array_1d<double,3>& cad_point_coords_in_deformed_configuration) override
    {
        Point<3> cad_point_in_deformed_configuration;
        mrAffectedPatch.EvaluateSurfacePoint(mParmeterValues, cad_point_in_deformed_configuration);
        
        cad_point_coords_in_deformed_configuration[0] = cad_point_in_deformed_configuration(0);
        cad_point_coords_in_deformed_configuration[1] = cad_point_in_deformed_configuration(1);
        cad_point_coords_in_deformed_configuration[2] = cad_point_in_deformed_configuration(2);        
    }
     
    // --------------------------------------------------------------------------
    Patch& GetAffectedPatch() override
    {
        return mrAffectedPatch;
    }  

    // --------------------------------------------------------------------------
    bool IsProjectedCADPointInsideVisiblePatchRegion() override
    {
        return mrAffectedPatch.IsPointInside(mParmeterValues);
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
        return "DistanceMinimizationCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DistanceMinimizationCondition";
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

  ///@}
  ///@name Private Operators
  ///@{

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    NodeType& mrFemNode;
    Patch& mrAffectedPatch;
    array_1d<double,2> mParmeterValues;
    array_1d<int,2> mParameterSpans;

    // ==============================================================================
    // Additional member variables
    // ==============================================================================
    std::vector<double> mNurbsFunctionValues;
    std::vector<ControlPoint*> mAffectedControlPoints;
    std::vector<int> mEquationIdsOfAffectedControlPoints;
    int mNumberOfLocalEquationIds;

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
  //      DistanceMinimizationCondition& operator=(DistanceMinimizationCondition const& rOther);

  /// Copy constructor.
  //      DistanceMinimizationCondition(DistanceMinimizationCondition const& rOther);

  ///@}

}; // Class DistanceMinimizationCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONDITION_DISTANCE_MINIMIZATION_H
