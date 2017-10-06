// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONDITION_DISPLACEMENT_MAPPING_H
#define RECONSTRUCTION_CONDITION_DISPLACEMENT_MAPPING_H

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

class DisplacementMappingCondition : public ReconstructionCondition
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element::GeometryType::IntegrationMethod IntegrationMethodType;
    typedef Element::GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /// Pointer definition of DisplacementMappingCondition
    KRATOS_CLASS_POINTER_DEFINITION(DisplacementMappingCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DisplacementMappingCondition( Element::GeometryType&  geometry,
                                  IntegrationMethodType int_method,
                                  int int_point_number,
                                  Patch& patch,
                                  array_1d<double,2> param_values,
                                  array_1d<double,2> param_spans )
    : mrGeometryContainingThisCondition( geometry ),
      mFemIntegrationMethod( int_method ),
      mIntegrationPointNumber( int_point_number ),
      mrAffectedPatch( patch ),
      mParmeterValues( param_values ),
      mParmeterSpans( param_spans )
    {
    }

    /// Destructor.
    virtual ~DisplacementMappingCondition()
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
        mrAffectedPatch.FlagAffectedControlPointsForReconstruction(  mParmeterSpans, mParmeterValues );
    }    
    
    // --------------------------------------------------------------------------
    void Initialize()
    {
        mIntegrationWeight = mrGeometryContainingThisCondition.IntegrationPoints(mFemIntegrationMethod)[mIntegrationPointNumber].Weight(); 
        
        mNurbsFunctionValues = mrAffectedPatch.EvaluateNURBSFunctions( mParmeterSpans, mParmeterValues );

        const Matrix& fem_shape_function_container = mrGeometryContainingThisCondition.ShapeFunctionsValues(mFemIntegrationMethod);
        mFEMFunctionValues = row( fem_shape_function_container, mIntegrationPointNumber);
        
        mAffectedControlPoints = mrAffectedPatch.GetPointersToAffectedControlPoints( mParmeterSpans, mParmeterValues );
        mEquationIdsOfAffectedControlPoints = mrAffectedPatch.GetEquationIdsOfAffectedControlPoints( mParmeterSpans, mParmeterValues );
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

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mIntegrationWeight * R_row * R_collumn;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mIntegrationWeight * R_row * R_collumn;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mIntegrationWeight * R_row * R_collumn;
            }
        }
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        int n_affected_fem_nodes =  mrGeometryContainingThisCondition.size();
        
        // Prepare vector of node displacements
        Vector node_displacements = ZeroVector(3*n_affected_fem_nodes);
        for(int itr  = 0; itr<n_affected_fem_nodes; itr++)
        {
            Vector node_disp = mrGeometryContainingThisCondition[itr].FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE);
            node_displacements[3*itr+0] = node_disp(0);
            node_displacements[3*itr+1] = node_disp(1);
            node_displacements[3*itr+2] = node_disp(2);
        }

        // Compute RHS
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIds; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPoints[row_itr];
            double R_row = mNurbsFunctionValues[row_itr];

            Vector corresponding_control_point_displacement = ZeroVector(3);
            corresponding_control_point_displacement[0] = mAffectedControlPoints[row_itr]->GetdX();
            corresponding_control_point_displacement[1] = mAffectedControlPoints[row_itr]->GetdY();
            corresponding_control_point_displacement[2] = mAffectedControlPoints[row_itr]->GetdZ();

            // Computation of RR*\hat{u_C}
            Vector rhs_contribution = ZeroVector(3);
            rhs_contribution(0) += R_row * R_row * corresponding_control_point_displacement[0];
            rhs_contribution(1) += R_row * R_row * corresponding_control_point_displacement[1];
            rhs_contribution(2) += R_row * R_row * corresponding_control_point_displacement[2]; 

            // Computation of -RN*\hat{u_F}
            for(int collumn_itr=0; collumn_itr<n_affected_fem_nodes; collumn_itr++)
            {
                double N_collumn = mFEMFunctionValues[collumn_itr];
                
                rhs_contribution(0) -= R_row * N_collumn * node_displacements[3*collumn_itr+0];
                rhs_contribution(1) -= R_row * N_collumn * node_displacements[3*collumn_itr+1];
                rhs_contribution(2) -= R_row * N_collumn * node_displacements[3*collumn_itr+2];
            }

            // Computation of complete RHS contribution: rhs_contribution = -integration_weight*( R*\hat{u_C} - N*\hat{u_F})R
            RHS( 3*row_id+0 ) -= mIntegrationWeight * rhs_contribution(0);
            RHS( 3*row_id+1 ) -= mIntegrationWeight * rhs_contribution(1);
            RHS( 3*row_id+2 ) -= mIntegrationWeight * rhs_contribution(2);             
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
        return "DisplacementMappingCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DisplacementMappingCondition";
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
    Element::GeometryType& mrGeometryContainingThisCondition;
    IntegrationMethodType mFemIntegrationMethod;
    int mIntegrationPointNumber;
    Patch& mrAffectedPatch;
    array_1d<double,2> mParmeterValues;
    array_1d<double,2>mParmeterSpans;

    // ==============================================================================
    // Additional member variables
    // ==============================================================================
    double mIntegrationWeight;
    std::vector<double> mNurbsFunctionValues;
    Vector mFEMFunctionValues; 
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
  //      DisplacementMappingCondition& operator=(DisplacementMappingCondition const& rOther);

  /// Copy constructor.
  //      DisplacementMappingCondition(DisplacementMappingCondition const& rOther);

  ///@}

}; // Class DisplacementMappingCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONDITION_DISPLACEMENT_MAPPING_H
