// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONSTRAINT_ZERO_DISPLACEMENT_H
#define RECONSTRUCTION_CONSTRAINT_ZERO_DISPLACEMENT_H

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

class ZeroDisplacementConstraint : public ReconstructionConstraint
{
public:
    ///@name Type Definitions
    ///@{   

    /// Pointer definition of ZeroDisplacementConstraint
    KRATOS_CLASS_POINTER_DEFINITION(ZeroDisplacementConstraint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ZeroDisplacementConstraint( BREPGaussPoint& coupling_gauss_point,
                                Patch& master_patch,
                                double penalty_factor )
    : mrCouplingGaussPoint( coupling_gauss_point ),
      mrAffectedMasterPatch( master_patch ),
      mPenaltyFactor( penalty_factor )
    {
        mLocationOnMaster = mrCouplingGaussPoint.GetLocationOnMasterInParameterSpace();
        mSpanOnMaster = mrAffectedMasterPatch.ComputeSurfaceKnotSpans( mLocationOnMaster );
    }

    /// Destructor.
    virtual ~ZeroDisplacementConstraint()
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
        mrAffectedMasterPatch.FlagAffectedControlPointsForReconstruction(  mSpanOnMaster, mLocationOnMaster );
    }

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        mIntegrationWeight = mrCouplingGaussPoint.GetWeight();
        
        mNurbsFunctionValuesOnMaster = mrAffectedMasterPatch.EvaluateNURBSFunctions( mSpanOnMaster, mLocationOnMaster );
        mEquationIdsOfAffectedControlPointsOnMaster = mrAffectedMasterPatch.GetEquationIdsOfAffectedControlPoints( mSpanOnMaster, mLocationOnMaster );
        mNumberOfLocalEquationIdsOnMaster = mEquationIdsOfAffectedControlPointsOnMaster.size();         
        
        // Compute Jacobian mJ1
        array_1d<double,2> tangent_on_master = mrCouplingGaussPoint.GetTangentOnMasterInParameterSpace();
        Matrix g_master = mrAffectedMasterPatch.ComputeBaseVectors( mSpanOnMaster, mLocationOnMaster );
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
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        array_1d<double,3> master_surface_displacement;
        mrAffectedMasterPatch.EvaluateSurfaceDisplacement( mLocationOnMaster, master_surface_displacement );
        
        for(int equation_itr=0; equation_itr<mNumberOfLocalEquationIdsOnMaster; equation_itr++)
        {
            int equation_id_master = mEquationIdsOfAffectedControlPointsOnMaster[equation_itr];
            double R_master = mNurbsFunctionValuesOnMaster[equation_itr];
            
            RHS[3*equation_id_master+0] -= mPenaltyFactor * mIntegrationWeight * mJ1 * ( master_surface_displacement[0] - 0.0 ) * R_master;
            RHS[3*equation_id_master+1] -= mPenaltyFactor * mIntegrationWeight * mJ1 * ( master_surface_displacement[1] - 0.0 ) * R_master;
            RHS[3*equation_id_master+2] -= mPenaltyFactor * mIntegrationWeight * mJ1 * ( master_surface_displacement[2] - 0.0 ) * R_master;           
        }
    }

    // --------------------------------------------------------------------------
    void Set( std::string identifier, double value ) override
    {
        if(identifier.compare("PENALTY_MULTIPLIER") == 0)
            mPenaltyFactor *= value;
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
        return "ZeroDisplacementConstraint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "ZeroDisplacementConstraint";
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
    double mPenaltyFactor = 0;

    // ==============================================================================
    // Additional member variables
    // ==============================================================================
    array_1d<double,2> mLocationOnMaster;
    array_1d<int,2> mSpanOnMaster;
    double mIntegrationWeight;
    std::vector<double> mNurbsFunctionValuesOnMaster;
    std::vector<int> mEquationIdsOfAffectedControlPointsOnMaster;
    int mNumberOfLocalEquationIdsOnMaster;
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
    //      ZeroDisplacementConstraint& operator=(ZeroDisplacementConstraint const& rOther);

    /// Copy constructor.
    //      ZeroDisplacementConstraint(ZeroDisplacementConstraint const& rOther);

    ///@}

}; // Class ZeroDisplacementConstraint

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONSTRAINT_ZERO_DISPLACEMENT_H
