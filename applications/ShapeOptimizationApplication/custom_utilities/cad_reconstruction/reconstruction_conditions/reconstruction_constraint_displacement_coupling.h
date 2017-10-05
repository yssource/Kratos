// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONSTRAINT_DISPLACEMENT_COUPLING_H
#define RECONSTRUCTION_CONSTRAINT_DISPLACEMENT_COUPLING_H

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
#include "reconstruction_constraint_base.h"

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

class DisplacementCouplingConstraint : public ReconstructionConstraint
{
public:
    ///@name Type Definitions
    ///@{   

    /// Pointer definition of DisplacementCouplingConstraint
    KRATOS_CLASS_POINTER_DEFINITION(DisplacementCouplingConstraint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DisplacementCouplingConstraint( BREPGaussPoint& coupling_gauss_point,
                                    Patch& master_patch,
                                    Patch& slave_patch,
                                    double penalty_factor )
    : mrCouplingGaussPoint( coupling_gauss_point ),
      mrAffectedMasterPatch( master_patch ),
      mrAffectedSlavePatch( slave_patch ),
      mPenaltyFactor( penalty_factor )
    {
        mLocationOnMaster = mrCouplingGaussPoint.GetLocationOnMasterInParameterSpace();
        mLocationOnSlave = mrCouplingGaussPoint.GetLocationOnSlaveInParameterSpace();        
    }

    /// Destructor.
    virtual ~DisplacementCouplingConstraint()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void FlagControlPointsRelevantForReconstruction() override
    {
        array_1d<int,2> unknown_parameter_span;
        unknown_parameter_span[0] = -1;
        unknown_parameter_span[1] = -1;

        mrAffectedMasterPatch.FlagAffectedControlPointsForReconstruction(  unknown_parameter_span, mLocationOnMaster );
        mrAffectedSlavePatch.FlagAffectedControlPointsForReconstruction(  unknown_parameter_span, mLocationOnSlave );      
    }

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        array_1d<double,2> tangent_on_master = mrCouplingGaussPoint.GetTangentOnMasterInParameterSpace();
        array_1d<int,2> unknown_parameter_span;
        unknown_parameter_span[0] = -1;
        unknown_parameter_span[1] = -1;
        
        mIntegrationWeight = mrCouplingGaussPoint.GetWeight();

        mNurbsFunctionValuesOnMaster = mrAffectedMasterPatch.EvaluateNURBSFunctions( unknown_parameter_span, mLocationOnMaster );
        mEquationIdsOfAffectedControlPointsOnMaster = mrAffectedMasterPatch.GetEquationIdsOfAffectedControlPoints( unknown_parameter_span, mLocationOnMaster );
        mNumberOfLocalEquationIdsOnMaster = mEquationIdsOfAffectedControlPointsOnMaster.size();  

        mNurbsFunctionValuesOnSlave = mrAffectedSlavePatch.EvaluateNURBSFunctions( unknown_parameter_span, mLocationOnSlave );
        mEquationIdsOfAffectedControlPointsOnSlave = mrAffectedSlavePatch.GetEquationIdsOfAffectedControlPoints( unknown_parameter_span, mLocationOnSlave );
        mNumberOfLocalEquationIdsOnSlave = mEquationIdsOfAffectedControlPointsOnSlave.size();          
        
        // Compute Jacobian mJ1
        Matrix g_master = mrAffectedMasterPatch.ComputeBaseVectors( unknown_parameter_span, mLocationOnMaster );
        Vector g1 = ZeroVector(3);
        g1(0) = g_master(0,0);
        g1(1) = g_master(1,0);
        g1(2) = g_master(2,0);
        Vector g2 = ZeroVector(3);
        g2(0) = g_master(0,1);
        g2(1) = g_master(1,1);
        g2(2) = g_master(2,1);
        mJ1 = norm_2( g1* tangent_on_master(0) + g2* tangent_on_master(1) );        
    }

    // --------------------------------------------------------------------------    
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) override
    {	
 	    // First we consider the relation Master-Master ( MM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];
            double R_row = mNurbsFunctionValuesOnMaster[row_itr];

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnMaster; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnMaster[collumn_itr];
                double R_collumn = mNurbsFunctionValuesOnMaster[collumn_itr];

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
            }
        }

 	    // Then we consider the relation Slave-Slave ( SS )
         for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
         {
             int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];
             double R_row = mNurbsFunctionValuesOnSlave[row_itr];
 
             for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
             {                
                 int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];
                 double R_collumn = mNurbsFunctionValuesOnSlave[collumn_itr];
 
                 LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                 LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                 LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
             }
         }        		

        // Then we consider the Master-Slave relation ( MS & SM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];
            double R_row = mNurbsFunctionValuesOnMaster[row_itr];

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];
                double R_collumn = mNurbsFunctionValuesOnSlave[collumn_itr];

                // MS
                LHS( 3*row_id+0, 3*collumn_id+0 ) -= mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                LHS( 3*row_id+1, 3*collumn_id+1 ) -= mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                LHS( 3*row_id+2, 3*collumn_id+2 ) -= mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                
                // SM
                LHS( 3*collumn_id+0, 3*row_id+0 ) -= mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                LHS( 3*collumn_id+1, 3*row_id+1 ) -= mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;
                LHS( 3*collumn_id+2, 3*row_id+2 ) -= mPenaltyFactor * mIntegrationWeight * mJ1 * R_row * R_collumn;                
            }
        }    
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        // array_1d<double,3> master_surface_displacement;
        // array_1d<double,3> slave_surface_displacement;

        // double patch_i.EvaluateSurfaceDisplacement( point_in_parameter_space, master_surface_displacement );
        // double patch_i.EvaluateSurfaceDisplacement( point_in_parameter_space, slave_surface_displacement );


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
        return "DisplacementCouplingConstraint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DisplacementCouplingConstraint";
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
    BREPGaussPoint& mrCouplingGaussPoint;
    Patch& mrAffectedMasterPatch;
    Patch& mrAffectedSlavePatch; 
    double mPenaltyFactor = 0;

    // ==============================================================================
    // Additional member variables
    // ==============================================================================
    double mIntegrationWeight;
    array_1d<double,2> mLocationOnMaster;
    array_1d<double,2> mLocationOnSlave;    
    std::vector<double> mNurbsFunctionValuesOnMaster;
    std::vector<int> mEquationIdsOfAffectedControlPointsOnMaster;
    int mNumberOfLocalEquationIdsOnMaster;
    std::vector<double> mNurbsFunctionValuesOnSlave;
    std::vector<int> mEquationIdsOfAffectedControlPointsOnSlave;
    int mNumberOfLocalEquationIdsOnSlave;
    double mJ1;     

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
    //      DisplacementCouplingConstraint& operator=(DisplacementCouplingConstraint const& rOther);

    /// Copy constructor.
    //      DisplacementCouplingConstraint(DisplacementCouplingConstraint const& rOther);

    ///@}

}; // Class DisplacementCouplingConstraint

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONSTRAINT_DISPLACEMENT_COUPLING_H
