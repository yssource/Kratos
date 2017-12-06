// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONSTRAINT_TANGENT_ENFORCEMENT_H
#define RECONSTRUCTION_CONSTRAINT_TANGENT_ENFORCEMENT_H

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

class TangentEnforcementConstraint : public ReconstructionConstraint
{
public:
    ///@name Type Definitions
    ///@{   

    /// Pointer definition of TangentEnforcementConstraint
    KRATOS_CLASS_POINTER_DEFINITION(TangentEnforcementConstraint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TangentEnforcementConstraint( BREPGaussPoint& coupling_gauss_point,
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
    virtual ~TangentEnforcementConstraint()
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
        Vector T1_der_m, T1_der_s, T2_der_m, T2_der_s, T3_der_m, T3_der_s;
        std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;
        std::vector<Vector> t1_der_r_m, t1_der_r_s, t2_der_r_m, t2_der_r_s, t3_der_r_m, t3_der_r_s;					
        std::vector<std::vector<Vector>> t1rs_m, t1rs_s, t2rs_m, t2rs_s, t3rs_m, t3rs_s;
        std::vector<std::vector<Vector>> t1_der_rs_m, t1_der_rs_s, t2_der_rs_m, t2_der_rs_s, t3_der_rs_m, t3_der_rs_s;

        // Compute geometric quantities
        mrAffectedMasterPatch.ComputeSecondVariationOfLocalCSY( mSpanOnMaster, 
                                                                mLocationOnMaster, 
                                                                mTangentOnMaster, 
                                                                T1_m, T2_m, T3_m, 
                                                                T1_der_m, T2_der_m, T3_der_m,
                                                                t1r_m, t2r_m, t3r_m,
                                                                t1_der_r_m, t2_der_r_m, t3_der_r_m,
                                                                t1rs_m, t2rs_m, t3rs_m,
                                                                t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
        mrAffectedSlavePatch.ComputeSecondVariationOfLocalCSY( mSpanOnSlave, 
                                                               mLocationOnSlave, 
                                                               mTangentOnSlave, 
                                                               T1_s, T2_s, T3_s, 
                                                               T1_der_s, T2_der_s, T3_der_s,
                                                               t1r_s, t2r_s, t3r_s,
                                                               t1_der_r_s, t2_der_r_s, t3_der_r_s,
                                                               t1rs_s, t2rs_s, t3rs_s,
                                                               t1_der_rs_s, t2_der_rs_s, t3_der_rs_s );

        double fac = ComputePreFactor( T3_m, T1_s );


        // First we consider the relation Master-Master ( MM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnMaster; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnMaster[collumn_itr];

                double term_1_x = inner_prod(t3r_m[3*collumn_itr+0],T1_s) * inner_prod(t3r_m[3*row_itr+0],T1_s);
                double term_1_y = inner_prod(t3r_m[3*collumn_itr+1],T1_s) * inner_prod(t3r_m[3*row_itr+1],T1_s);
                double term_1_z = inner_prod(t3r_m[3*collumn_itr+2],T1_s) * inner_prod(t3r_m[3*row_itr+2],T1_s);

                double term_2_x = fac * inner_prod(t3rs_m[3*row_itr+0][3*collumn_itr+0],T1_s);
                double term_2_y = fac * inner_prod(t3rs_m[3*row_itr+1][3*collumn_itr+1],T1_s);
                double term_2_z = fac * inner_prod(t3rs_m[3*row_itr+2][3*collumn_itr+2],T1_s);

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_x + term_2_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_y + term_2_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_z + term_2_z );	
            }
        }

        // Then we consider the relation Slave-Slave ( SS )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr]; 

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                double term_1_x = inner_prod(t1r_s[3*collumn_itr+0],T3_m) * inner_prod(T3_m, t1r_s[3*row_itr+0]);
                double term_1_y = inner_prod(t1r_s[3*collumn_itr+1],T3_m) * inner_prod(T3_m, t1r_s[3*row_itr+1]);
                double term_1_z = inner_prod(t1r_s[3*collumn_itr+2],T3_m) * inner_prod(T3_m, t1r_s[3*row_itr+2]);

                double term_2_x = fac * inner_prod(T3_m, t1rs_s[3*row_itr+0][3*collumn_itr+0]);
                double term_2_y = fac * inner_prod(T3_m, t1rs_s[3*row_itr+1][3*collumn_itr+1]);
                double term_2_z = fac * inner_prod(T3_m, t1rs_s[3*row_itr+2][3*collumn_itr+2]);

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_x + term_2_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_y + term_2_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_z + term_2_z );	
            }
        }        

        // Then we consider the Master-Slave relation ( MS )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];
              
                double term_1_x = inner_prod(T3_m,t1r_s[3*collumn_itr+0]) * inner_prod(t3r_m[3*row_itr+0],T1_s);
                double term_1_y = inner_prod(T3_m,t1r_s[3*collumn_itr+1]) * inner_prod(t3r_m[3*row_itr+1],T1_s);
                double term_1_z = inner_prod(T3_m,t1r_s[3*collumn_itr+2]) * inner_prod(t3r_m[3*row_itr+2],T1_s);

                double term_2_x = fac * inner_prod(t3r_m[3*row_itr+0],t1r_s[3*collumn_itr+0]);
                double term_2_y = fac * inner_prod(t3r_m[3*row_itr+1],t1r_s[3*collumn_itr+1]);
                double term_2_z = fac * inner_prod(t3r_m[3*row_itr+2],t1r_s[3*collumn_itr+2]);

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_x + term_2_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_y + term_2_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_z + term_2_z );
            }
        }  

        // Then we consider the Slave-Master relation ( SM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnMaster; collumn_itr++)
            {                
                int collumn_id = mEquationIdsOfAffectedControlPointsOnMaster[collumn_itr];
              
                double term_1_y = inner_prod(t3r_m[3*collumn_itr+1], T1_s) * inner_prod(T3_m, t1r_s[3*row_itr+1]);
                double term_1_z = inner_prod(t3r_m[3*collumn_itr+2], T1_s) * inner_prod(T3_m, t1r_s[3*row_itr+2]);
                double term_1_x = inner_prod(t3r_m[3*collumn_itr+0], T1_s) * inner_prod(T3_m, t1r_s[3*row_itr+0]);

                double term_2_x = fac * inner_prod(t3r_m[3*collumn_itr+0], t1r_s[3*row_itr+0]);
                double term_2_y = fac * inner_prod(t3r_m[3*collumn_itr+1], t1r_s[3*row_itr+1]);
                double term_2_z = fac * inner_prod(t3r_m[3*collumn_itr+2], t1r_s[3*row_itr+2]);

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_x + term_2_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_y + term_2_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( term_1_z + term_2_z );	
            }
        }
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        // Variables needed later
        Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
        Vector T1_der_m, T1_der_s, T2_der_m, T2_der_s, T3_der_m, T3_der_s;
        std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;
        std::vector<Vector> t1_der_r_m, t1_der_r_s, t2_der_r_m, t2_der_r_s, t3_der_r_m, t3_der_r_s;					
        std::vector<std::vector<Vector>> t1rs_m, t1rs_s, t2rs_m, t2rs_s, t3rs_m, t3rs_s;
        std::vector<std::vector<Vector>> t1_der_rs_m, t1_der_rs_s, t2_der_rs_m, t2_der_rs_s, t3_der_rs_m, t3_der_rs_s;

        // Compute geometric quantities
        mrAffectedMasterPatch.ComputeSecondVariationOfLocalCSY( mSpanOnMaster, 
                                                                mLocationOnMaster, 
                                                                mTangentOnMaster, 
                                                                T1_m, T2_m, T3_m, 
                                                                T1_der_m, T2_der_m, T3_der_m,
                                                                t1r_m, t2r_m, t3r_m,
                                                                t1_der_r_m, t2_der_r_m, t3_der_r_m,
                                                                t1rs_m, t2rs_m, t3rs_m,
                                                                t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
        mrAffectedSlavePatch.ComputeSecondVariationOfLocalCSY( mSpanOnSlave, 
                                                               mLocationOnSlave, 
                                                               mTangentOnSlave, 
                                                               T1_s, T2_s, T3_s, 
                                                               T1_der_s, T2_der_s, T3_der_s,
                                                               t1r_s, t2r_s, t3r_s,
                                                               t1_der_r_s, t2_der_r_s, t3_der_r_s,
                                                               t1rs_s, t2rs_s, t3rs_s,
                                                               t1_der_rs_s, t2_der_rs_s, t3_der_rs_s );

        double fac = ComputePreFactor( T3_m, T1_s );


        // First we consider the relation Master-Master ( MM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            RHS( 3*row_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * fac * inner_prod(t3r_m[3*row_itr+0],T1_s);
            RHS( 3*row_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * fac * inner_prod(t3r_m[3*row_itr+1],T1_s);
            RHS( 3*row_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * fac * inner_prod(t3r_m[3*row_itr+2],T1_s);
        }

        // The we consider the relation Slave-Slave ( SS )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            RHS( 3*row_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * fac * inner_prod(t1r_s[3*row_itr+0],T3_m);
            RHS( 3*row_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * fac * inner_prod(t1r_s[3*row_itr+1],T3_m);
            RHS( 3*row_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * fac * inner_prod(t1r_s[3*row_itr+2],T3_m);            
        }
    }

    // --------------------------------------------------------------------------
    double ComputePreFactor( Vector T3_m, Vector T1_s )
    {
        return - inner_prod(T3_m, T1_s);
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
        return "TangentEnforcementConstraint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "TangentEnforcementConstraint";
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
    //      TangentEnforcementConstraint& operator=(TangentEnforcementConstraint const& rOther);

    /// Copy constructor.
    //      TangentEnforcementConstraint(TangentEnforcementConstraint const& rOther);

    ///@}

}; // Class TangentEnforcementConstraint

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONSTRAINT_TANGENT_ENFORCEMENT_H
