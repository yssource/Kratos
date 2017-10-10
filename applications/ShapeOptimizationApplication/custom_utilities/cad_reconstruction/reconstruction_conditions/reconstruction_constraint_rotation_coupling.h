// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONSTRAINT_ROTATION_COUPLING_H
#define RECONSTRUCTION_CONSTRAINT_ROTATION_COUPLING_H

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

class RotationCouplingConstraint : public ReconstructionConstraint
{
public:
    ///@name Type Definitions
    ///@{   

    /// Pointer definition of RotationCouplingConstraint
    KRATOS_CLASS_POINTER_DEFINITION(RotationCouplingConstraint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RotationCouplingConstraint( BREPGaussPoint& coupling_gauss_point,
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
        mSpanOnMaster = mrAffectedMasterPatch.ComputeSurfaceKnotSpans( mLocationOnMaster );
        mSpanOnSlave = mrAffectedSlavePatch.ComputeSurfaceKnotSpans( mLocationOnSlave );        
    }

    /// Destructor.
    virtual ~RotationCouplingConstraint()
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
        mrAffectedSlavePatch.FlagAffectedControlPointsForReconstruction(  mSpanOnSlave, mLocationOnSlave );      
    }

    // --------------------------------------------------------------------------
    void Initialize() override
    {        
        mIntegrationWeight = mrCouplingGaussPoint.GetWeight();
        
        mTangentOnMaster = mrCouplingGaussPoint.GetTangentOnMasterInParameterSpace();
        mTangentOnSlave = mrCouplingGaussPoint.GetTangentOnSlaveInParameterSpace();

        mEquationIdsOfAffectedControlPointsOnMaster = mrAffectedMasterPatch.GetEquationIdsOfAffectedControlPoints( mSpanOnMaster, mLocationOnMaster );
        mNumberOfLocalEquationIdsOnMaster = mEquationIdsOfAffectedControlPointsOnMaster.size();  
        
        mEquationIdsOfAffectedControlPointsOnSlave = mrAffectedSlavePatch.GetEquationIdsOfAffectedControlPoints( mSpanOnSlave, mLocationOnSlave );
        mNumberOfLocalEquationIdsOnSlave = mEquationIdsOfAffectedControlPointsOnSlave.size();          
        
        // Compute Jacobian mJ1
        Matrix g_master = mrAffectedMasterPatch.ComputeBaseVectors( mSpanOnMaster, mLocationOnMaster );
        Vector g1 = ZeroVector(3);
        g1(0) = g_master(0,0);
        g1(1) = g_master(1,0);
        g1(2) = g_master(2,0);
        Vector g2 = ZeroVector(3);
        g2(0) = g_master(0,1);
        g2(1) = g_master(1,1);
        g2(2) = g_master(2,1);
        mJ1 = norm_2( g1* mTangentOnMaster(0) + g2* mTangentOnMaster(1) );        
    }

    // --------------------------------------------------------------------------    
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) override
    {
        // Variables needed later
        Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
        std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;				

        // Compute geometric quantities
        mrAffectedMasterPatch.ComputeVariationOfLocalCSY( mSpanOnMaster, mLocationOnMaster, mTangentOnMaster, T1_m, T2_m, T3_m, t1r_m, t2r_m, t3r_m );
        mrAffectedSlavePatch.ComputeVariationOfLocalCSY( mSpanOnSlave, mLocationOnSlave, mTangentOnSlave, T1_s, T2_s, T3_s, t1r_s, t2r_s, t3r_s );

        // Check if master and slave tangent point in same direction. If yes, we have to subtract in the following.
        int sign_factor = 1;
        if( inner_prod(T2_m,T2_s) > 0 )
            sign_factor = -1;

        // First we consider the relation Master-Master ( MM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_mx_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_my_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_mz_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_T2_mx_row = inner_prod(omega_mx_row,T2_m);
            double omega_T2_my_row = inner_prod(omega_my_row,T2_m);
            double omega_T2_mz_row = inner_prod(omega_mz_row,T2_m);            

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnMaster; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnMaster[collumn_itr];

                Vector omega_mx_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+0]);
                Vector omega_my_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+1]);
                Vector omega_mz_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+2]);
                double omega_T2_mx_coll = inner_prod(omega_mx_coll,T2_m);
                double omega_T2_my_coll = inner_prod(omega_my_coll,T2_m);
                double omega_T2_mz_coll = inner_prod(omega_mz_coll,T2_m);

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_mx_row * omega_T2_mx_coll;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_my_row * omega_T2_my_coll;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_mz_row * omega_T2_mz_coll;
            }
        }

        // Then we consider the relation Slave-Slave ( SS )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            Vector omega_sx_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+0]);
            Vector omega_sy_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+1]);
            Vector omega_sz_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+2]);
            double omega_T2_sx_row = inner_prod(omega_sx_row,T2_s);
            double omega_T2_sy_row = inner_prod(omega_sy_row,T2_s);
            double omega_T2_sz_row = inner_prod(omega_sz_row,T2_s);            

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                Vector omega_sx_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+0]);
                Vector omega_sy_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+1]);
                Vector omega_sz_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+2]);
                double omega_T2_sx_coll = inner_prod(omega_sx_coll,T2_s);
                double omega_T2_sy_coll = inner_prod(omega_sy_coll,T2_s);
                double omega_T2_sz_coll = inner_prod(omega_sz_coll,T2_s);

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_sx_row * omega_T2_sx_coll;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_sy_row * omega_T2_sy_coll;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_sz_row * omega_T2_sz_coll;
            }
        } 
		
        // Then we consider the Master-Slave relation ( MS & SM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_mx = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_my = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_mz = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_T2_mx = inner_prod(omega_mx,T2_m);
            double omega_T2_my = inner_prod(omega_my,T2_m);
            double omega_T2_mz = inner_prod(omega_mz,T2_m);

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                Vector omega_sx = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+0]);
                Vector omega_sy = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+1]);
                Vector omega_sz = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+2]);
                double omega_T2_sx = inner_prod(omega_sx,T2_s);
                double omega_T2_sy = inner_prod(omega_sy,T2_s);
                double omega_T2_sz = inner_prod(omega_sz,T2_s);

                // MS
                LHS( 3*row_id+0, 3*collumn_id+0 ) += sign_factor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_mx * omega_T2_sx;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += sign_factor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_my * omega_T2_sy;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += sign_factor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_mz * omega_T2_sz;
                
                // SM
                LHS( 3*collumn_id+0, 3*row_id+0 ) += sign_factor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_mx * omega_T2_sx;
                LHS( 3*collumn_id+1, 3*row_id+1 ) += sign_factor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_my * omega_T2_sy;
                LHS( 3*collumn_id+2, 3*row_id+2 ) += sign_factor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_mz * omega_T2_sz;                
            }
        }  
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        // To be implemented
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
        return "RotationCouplingConstraint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "RotationCouplingConstraint";
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
    array_1d<double,2> mLocationOnMaster;
    array_1d<double,2> mLocationOnSlave;
    array_1d<int,2> mSpanOnMaster;
    array_1d<int,2> mSpanOnSlave;    
    double mIntegrationWeight;
    array_1d<double,2> mTangentOnMaster;
    array_1d<double,2> mTangentOnSlave;
    std::vector<int> mEquationIdsOfAffectedControlPointsOnMaster;
    int mNumberOfLocalEquationIdsOnMaster;
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
    //      RotationCouplingConstraint& operator=(RotationCouplingConstraint const& rOther);

    /// Copy constructor.
    //      RotationCouplingConstraint(RotationCouplingConstraint const& rOther);

    ///@}

}; // Class RotationCouplingConstraint

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONSTRAINT_ROTATION_COUPLING_H
