// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONSTRAINT_ROTATION_TARGET_H
#define RECONSTRUCTION_CONSTRAINT_ROTATION_TARGET_H

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

class RotationTargetConstraint : public ReconstructionConstraint
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RotationTargetConstraint
    KRATOS_CLASS_POINTER_DEFINITION(RotationTargetConstraint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RotationTargetConstraint( BREPGaussPoint& coupling_gauss_point,
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
    virtual ~RotationTargetConstraint()
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
        Vector g1_master = ZeroVector(3);
        g1_master(0) = g_master(0,0);
        g1_master(1) = g_master(1,0);
        g1_master(2) = g_master(2,0);
        Vector g2_master = ZeroVector(3);
        g2_master(0) = g_master(0,1);
        g2_master(1) = g_master(1,1);
        g2_master(2) = g_master(2,1);
        mJ1 = norm_2( g1_master* mTangentOnMaster(0) + g2_master* mTangentOnMaster(1) );

        // Compute initial coordinate system
        mrAffectedMasterPatch.ComputeLocalCSY( mSpanOnMaster, mLocationOnMaster, mTangentOnMaster, T1_m, T2_m, T3_m );
        mrAffectedSlavePatch.ComputeLocalCSY( mSpanOnSlave, mLocationOnSlave, mTangentOnSlave, T1_s, T2_s, T3_s );

        // Check if master and slave tangent point in same direction. If yes, we have to subtract in the following.
        mSignFactor = 1;
        if( inner_prod(T2_m,T2_s) > 0 )
            mSignFactor = -1;

        KRATOS_WATCH(mSignFactor)

        // Compute target rotation
        Vector w = mSignFactor*T3_s - T3_m;
        Vector omega = MathUtils<double>::CrossProduct(T3_m,w);
        mTargetRotation = std::asin( inner_prod(omega, T2_m) );
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) override
    {
        // Variables needed later
        Vector t1_m, t1_s, t2_m, t2_s, t3_m, t3_s;
        std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;

        // Compute geometric quantities
        mrAffectedMasterPatch.ComputeVariationOfLocalCSY( mSpanOnMaster, mLocationOnMaster, mTangentOnMaster, t1_m, t2_m, t3_m, t1r_m, t2r_m, t3r_m );
        mrAffectedSlavePatch.ComputeVariationOfLocalCSY( mSpanOnSlave, mLocationOnSlave, mTangentOnSlave, t1_s, t2_s, t3_s, t1r_s, t2r_s, t3r_s );

        // Compute difference between rotation of master and slave
        Vector w_m = t3_m - T3_m;
        Vector w_s = t3_s - T3_s;
        Vector omega_m = MathUtils<double>::CrossProduct(T3_m,w_m);
        Vector omega_s = MathUtils<double>::CrossProduct(T3_s,w_s);
        double tmp_inner_prod_m = inner_prod(omega_m, T2_m);
        double tmp_inner_prod_s = inner_prod(omega_s, T2_s);

        // First we consider the relation Master-Master ( MM )
        double divisor_m = std::sqrt((1-tmp_inner_prod_m*tmp_inner_prod_m));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_r_T2_x = inner_prod(omega_r_x,T2_m) / divisor_m;
            double omega_r_T2_y = inner_prod(omega_r_y,T2_m) / divisor_m;
            double omega_r_T2_z = inner_prod(omega_r_z,T2_m) / divisor_m;

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnMaster; collumn_itr++)
            {
                int collumn_id = mEquationIdsOfAffectedControlPointsOnMaster[collumn_itr];

                Vector omega_s_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+0]);
                Vector omega_s_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+1]);
                Vector omega_s_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+2]);
                double omega_s_T2_x = inner_prod(omega_s_x,T2_m) / divisor_m;
                double omega_s_T2_y = inner_prod(omega_s_y,T2_m) / divisor_m;
                double omega_s_T2_z = inner_prod(omega_s_z,T2_m) / divisor_m;

                double second_order_term_x = 0.0; //(current_rotation_difference - mTargetRotation) * ...;
                double second_order_term_y = 0.0; //(current_rotation_difference - mTargetRotation) * ...;
                double second_order_term_z = 0.0; //(current_rotation_difference - mTargetRotation) * ...;

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_r_T2_x * omega_s_T2_x +  second_order_term_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_r_T2_y * omega_s_T2_y +  second_order_term_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_r_T2_z * omega_s_T2_z +  second_order_term_z );
            }
        }

        // Then we consider the relation Slave-Slave ( SS )
        double divisor_s = std::sqrt((1-tmp_inner_prod_s*tmp_inner_prod_s));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+2]);
            double omega_r_T2_x = inner_prod(omega_r_x,T2_s) / divisor_s;
            double omega_r_T2_y = inner_prod(omega_r_y,T2_s) / divisor_s;
            double omega_r_T2_z = inner_prod(omega_r_z,T2_s) / divisor_s;

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                Vector omega_s_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+0]);
                Vector omega_s_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+1]);
                Vector omega_s_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+2]);
                double omega_s_T2_x = inner_prod(omega_s_x,T2_s) / divisor_s;
                double omega_s_T2_y = inner_prod(omega_s_y,T2_s) / divisor_s;
                double omega_s_T2_z = inner_prod(omega_s_z,T2_s) / divisor_s;

                double second_order_term_x = 0.0; //(current_rotation_difference - mTargetRotation) ...;
                double second_order_term_y = 0.0; //(current_rotation_difference - mTargetRotation) ...;
                double second_order_term_z = 0.0; //(current_rotation_difference - mTargetRotation) ...;

                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_r_T2_x * omega_s_T2_x +  second_order_term_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_r_T2_y * omega_s_T2_y +  second_order_term_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_r_T2_z * omega_s_T2_z +  second_order_term_z );
            }
        }

        // Then we consider the Master-Slave relation ( MS & SM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_r_T2_x = inner_prod(omega_r_x,T2_m) / divisor_m;
            double omega_r_T2_y = inner_prod(omega_r_y,T2_m) / divisor_m;
            double omega_r_T2_z = inner_prod(omega_r_z,T2_m) / divisor_m;

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                Vector omega_s_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+0]);
                Vector omega_s_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+1]);
                Vector omega_s_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+2]);
                double omega_s_T2_x = inner_prod(omega_s_x,T2_s) / divisor_s;
                double omega_s_T2_y = inner_prod(omega_s_y,T2_s) / divisor_s;
                double omega_s_T2_z = inner_prod(omega_s_z,T2_s) / divisor_s;

                // MS
                LHS( 3*row_id+0, 3*collumn_id+0 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_r_T2_x * omega_s_T2_x;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_r_T2_y * omega_s_T2_y;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_r_T2_z * omega_s_T2_z;

                // SM
                LHS( 3*collumn_id+0, 3*row_id+0 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_r_T2_x * omega_s_T2_x;
                LHS( 3*collumn_id+1, 3*row_id+1 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_r_T2_y * omega_s_T2_y;
                LHS( 3*collumn_id+2, 3*row_id+2 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_r_T2_z * omega_s_T2_z;
            }
        }
    }

    // // Working version but at the very edges, there is problems (this approach anyways is mathimatically not very consistent)
    // // --------------------------------------------------------------------------
    // void ComputeAndAddRHSContribution( Vector& RHS ) override
    // {
    //     // Variables needed later
    //     Vector t1_m, t1_s, t2_m, t2_s, t3_m, t3_s;
    //     std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;

    //     // Compute geometric quantities
    //     mrAffectedMasterPatch.ComputeVariationOfLocalCSY( mSpanOnMaster, mLocationOnMaster, mTangentOnMaster, t1_m, t2_m, t3_m, t1r_m, t2r_m, t3r_m );
    //     mrAffectedSlavePatch.ComputeVariationOfLocalCSY( mSpanOnSlave, mLocationOnSlave, mTangentOnSlave, t1_s, t2_s, t3_s, t1r_s, t2r_s, t3r_s );

    //     // Check if master and slave tangent point in same direction. If yes, we have to subtract in the following.
    //     int mSignFactor = 1;
    //     if( inner_prod(T2_m,T2_s) > 0 )
    //         mSignFactor = -1;

    //     Vector w = mSignFactor*t3_s - t3_m;
    //     Vector omega = MathUtils<double>::CrossProduct(t3_m,w);
    //     double target_rotation = std::asin( inner_prod(omega, t2_m) );
    //     double target_rotation_in_deg = target_rotation*180/3.14;
    //     KRATOS_WATCH(target_rotation_in_deg)

    //     // Master contribution
    //     for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
    //     {
    //         int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

    //         Vector omega_mx_row = MathUtils<double>::CrossProduct(t3_m,t3r_m[3*row_itr+0]);
    //         Vector omega_my_row = MathUtils<double>::CrossProduct(t3_m,t3r_m[3*row_itr+1]);
    //         Vector omega_mz_row = MathUtils<double>::CrossProduct(t3_m,t3r_m[3*row_itr+2]);
    //         double omega_t2_mx_row = inner_prod(omega_mx_row,t2_m);
    //         double omega_t2_my_row = inner_prod(omega_my_row,t2_m);
    //         double omega_t2_mz_row = inner_prod(omega_mz_row,t2_m);

    //         RHS[3*row_id+0] -= mPenaltyFactor * mIntegrationWeight * mJ1 * ( - target_rotation) * omega_t2_mx_row;
    //         RHS[3*row_id+1] -= mPenaltyFactor * mIntegrationWeight * mJ1 * ( - target_rotation) * omega_t2_my_row;
    //         RHS[3*row_id+2] -= mPenaltyFactor * mIntegrationWeight * mJ1 * ( - target_rotation) * omega_t2_mz_row;
    //     }

    //     // Slave contribution
    //     for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
    //     {
    //         int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

    //         Vector omega_sx_row = MathUtils<double>::CrossProduct(t3_s,t3r_s[3*row_itr+0]);
    //         Vector omega_sy_row = MathUtils<double>::CrossProduct(t3_s,t3r_s[3*row_itr+1]);
    //         Vector omega_sz_row = MathUtils<double>::CrossProduct(t3_s,t3r_s[3*row_itr+2]);
    //         double omega_t2_sx_row = inner_prod(omega_sx_row,t2_s);
    //         double omega_t2_sy_row = inner_prod(omega_sy_row,t2_s);
    //         double omega_t2_sz_row = inner_prod(omega_sz_row,t2_s);

    //         RHS[3*row_id+0] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * ( - target_rotation ) * omega_t2_sx_row;
    //         RHS[3*row_id+1] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * ( - target_rotation ) * omega_t2_sy_row;
    //         RHS[3*row_id+2] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * ( - target_rotation ) * omega_t2_sz_row;
    //     }
    // }

    // Working version but after several iterations, surface gets bad in in-plane direction!
    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
        // Variables needed later
        Vector t1_m, t1_s, t2_m, t2_s, t3_m, t3_s;
        std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;

        // Compute geometric quantities
        mrAffectedMasterPatch.ComputeVariationOfLocalCSY( mSpanOnMaster, mLocationOnMaster, mTangentOnMaster, t1_m, t2_m, t3_m, t1r_m, t2r_m, t3r_m );
        mrAffectedSlavePatch.ComputeVariationOfLocalCSY( mSpanOnSlave, mLocationOnSlave, mTangentOnSlave, t1_s, t2_s, t3_s, t1r_s, t2r_s, t3r_s );

        // Compute difference between rotation of master and slave
        Vector w_m = t3_m - T3_m;
        Vector w_s = t3_s - T3_s;
        Vector omega_m = MathUtils<double>::CrossProduct(T3_m,w_m);
        Vector omega_s = MathUtils<double>::CrossProduct(T3_s,w_s);
        double tmp_inner_prod_m = inner_prod(omega_m, T2_m);
        double tmp_inner_prod_s = inner_prod(omega_s, T2_s);
        double current_rotation_m = std::asin( tmp_inner_prod_m );
        double current_rotation_s = std::asin( tmp_inner_prod_s );
        double current_rotation_difference = current_rotation_m + mSignFactor*current_rotation_s;

        double target_rotation_in_deg = mTargetRotation*180/3.14;
        double current_rotation_difference_in_deg = current_rotation_difference*180/3.14;

        KRATOS_WATCH("--------------")
        // KRATOS_WATCH(w_m)
        // KRATOS_WATCH(w_s)
        // KRATOS_WATCH(tmp_inner_prod_m)
        // KRATOS_WATCH(tmp_inner_prod_s)
        KRATOS_WATCH(target_rotation_in_deg)
        KRATOS_WATCH(current_rotation_difference_in_deg)

        // Master contribution
        double divisor_m = std::sqrt((1-tmp_inner_prod_m*tmp_inner_prod_m));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_r_T2_x = inner_prod(omega_r_x,T2_m) / divisor_m;
            double omega_r_T2_y = inner_prod(omega_r_y,T2_m) / divisor_m;
            double omega_r_T2_z = inner_prod(omega_r_z,T2_m) / divisor_m;

            RHS[3*row_id+0] -= mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_r_T2_x;
            RHS[3*row_id+1] -= mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_r_T2_y;
            RHS[3*row_id+2] -= mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_r_T2_z;
        }

        // Slave contribution
        double divisor_s = std::sqrt((1-tmp_inner_prod_s*tmp_inner_prod_s));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+2]);
            double omega_r_T2_x = inner_prod(omega_r_x,T2_s) / divisor_s;
            double omega_r_T2_y = inner_prod(omega_r_y,T2_s) / divisor_s;
            double omega_r_T2_z = inner_prod(omega_r_z,T2_s) / divisor_s;

            RHS[3*row_id+0] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_r_T2_x;
            RHS[3*row_id+1] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_r_T2_y;
            RHS[3*row_id+2] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_r_T2_z;
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
        return "RotationTargetConstraint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "RotationTargetConstraint";
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
    Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
    double mSignFactor;
    double mTargetRotation;

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
    //      RotationTargetConstraint& operator=(RotationTargetConstraint const& rOther);

    /// Copy constructor.
    //      RotationTargetConstraint(RotationTargetConstraint const& rOther);

    ///@}

}; // Class RotationTargetConstraint

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONSTRAINT_ROTATION_TARGET_H
