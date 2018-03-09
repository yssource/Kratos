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
        Vector t1_der_m, t1_der_s, t2_der_m, t2_der_s, t3_der_m, t3_der_s;
        std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;
        std::vector<Vector> t1_der_r_m, t1_der_r_s, t2_der_r_m, t2_der_r_s, t3_der_r_m, t3_der_r_s;
        std::vector<std::vector<Vector>> t1rs_m, t1rs_s, t2rs_m, t2rs_s, t3rs_m, t3rs_s;
        std::vector<std::vector<Vector>> t1_der_rs_m, t1_der_rs_s, t2_der_rs_m, t2_der_rs_s, t3_der_rs_m, t3_der_rs_s;

        // Compute geometric quantities
        mrAffectedMasterPatch.ComputeSecondVariationOfLocalCSY( mSpanOnMaster,
                                                                mLocationOnMaster,
                                                                mTangentOnMaster,
                                                                t1_m, t2_m, t3_m,
                                                                t1_der_m, t2_der_m, t3_der_m,
                                                                t1r_m, t2r_m, t3r_m,
                                                                t1_der_r_m, t2_der_r_m, t3_der_r_m,
                                                                t1rs_m, t2rs_m, t3rs_m,
                                                                t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
        mrAffectedSlavePatch.ComputeSecondVariationOfLocalCSY( mSpanOnSlave,
                                                               mLocationOnSlave,
                                                               mTangentOnSlave,
                                                               t1_s, t2_s, t3_s,
                                                               t1_der_s, t2_der_s, t3_der_s,
                                                               t1r_s, t2r_s, t3r_s,
                                                               t1_der_r_s, t2_der_r_s, t3_der_r_s,
                                                               t1rs_s, t2rs_s, t3rs_s,
                                                               t1_der_rs_s, t2_der_rs_s, t3_der_rs_s );


        // // Variables needed later
        // Vector t1_m, t1_s, t2_m, t2_s, t3_m, t3_s;
        // std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;

        // // Compute geometric quantities
        // mrAffectedMasterPatch.ComputeVariationOfLocalCSY( mSpanOnMaster, mLocationOnMaster, mTangentOnMaster, t1_m, t2_m, t3_m, t1r_m, t2r_m, t3r_m );
        // mrAffectedSlavePatch.ComputeVariationOfLocalCSY( mSpanOnSlave, mLocationOnSlave, mTangentOnSlave, t1_s, t2_s, t3_s, t1r_s, t2r_s, t3r_s );

        // Compute difference between rotation of master and slave
        Vector w_m = t3_m - T3_m;
        Vector w_s = t3_s - T3_s;
        Vector omega_m = MathUtils<double>::CrossProduct(T3_m,w_m);
        Vector omega_s = MathUtils<double>::CrossProduct(T3_s,w_s);
        double omega_dot_T2_m = inner_prod(omega_m, T2_m);
        double omega_dot_T2_s = inner_prod(omega_s, T2_s);
        double current_rotation_m = std::asin( omega_dot_T2_m );
        double current_rotation_s = std::asin( omega_dot_T2_s );
        double current_rotation_difference = current_rotation_m + mSignFactor*current_rotation_s;

        // First we consider the relation Master-Master ( MM )
        double divisor_m_1 = std::sqrt((1-omega_dot_T2_m*omega_dot_T2_m));
        double divisor_m_2 = std::pow((1-omega_dot_T2_m*omega_dot_T2_m),(3.0/2.0));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_r_x_dot_T2 = inner_prod(omega_r_x,T2_m);
            double omega_r_y_dot_T2 = inner_prod(omega_r_y,T2_m);
            double omega_r_z_dot_T2 = inner_prod(omega_r_z,T2_m);
            double omega_T2_r_x = omega_r_x_dot_T2 / divisor_m_1;
            double omega_T2_r_y = omega_r_y_dot_T2 / divisor_m_1;
            double omega_T2_r_z = omega_r_z_dot_T2 / divisor_m_1;

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnMaster; collumn_itr++)
            {
                int collumn_id = mEquationIdsOfAffectedControlPointsOnMaster[collumn_itr];

                Vector omega_s_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+0]);
                Vector omega_s_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+1]);
                Vector omega_s_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*collumn_itr+2]);
                double omega_s_x_dot_T2 = inner_prod(omega_s_x,T2_m);
                double omega_s_y_dot_T2 = inner_prod(omega_s_y,T2_m);
                double omega_s_z_dot_T2 = inner_prod(omega_s_z,T2_m);
                double omega_T2_s_x = omega_s_x_dot_T2 / divisor_m_1;
                double omega_T2_s_y = omega_s_y_dot_T2 / divisor_m_1;
                double omega_T2_s_z = omega_s_z_dot_T2 / divisor_m_1;

                // First term of second derivative part
                Vector omega_rs_x = MathUtils<double>::CrossProduct(T3_m,t3rs_m[3*row_itr+0][3*collumn_itr+0]);
                Vector omega_rs_y = MathUtils<double>::CrossProduct(T3_m,t3rs_m[3*row_itr+1][3*collumn_itr+1]);
                Vector omega_rs_z = MathUtils<double>::CrossProduct(T3_m,t3rs_m[3*row_itr+2][3*collumn_itr+2]);

                double omega_T2_rs_x = inner_prod(omega_rs_x,T2_m) / divisor_m_1;
                double omega_T2_rs_y = inner_prod(omega_rs_y,T2_m) / divisor_m_1;
                double omega_T2_rs_z = inner_prod(omega_rs_z,T2_m) / divisor_m_1;

                // Second term of second derivative part
                omega_T2_rs_x += (omega_r_x_dot_T2 * omega_dot_T2_m * omega_s_x_dot_T2) / divisor_m_2;
                omega_T2_rs_y += (omega_r_y_dot_T2 * omega_dot_T2_m * omega_s_y_dot_T2) / divisor_m_2;
                omega_T2_rs_z += (omega_r_z_dot_T2 * omega_dot_T2_m * omega_s_z_dot_T2) / divisor_m_2;

                // Complete second derivative part
                double second_order_term_x = omega_T2_rs_x * (current_rotation_difference - mTargetRotation);
                double second_order_term_y = omega_T2_rs_y * (current_rotation_difference - mTargetRotation);
                double second_order_term_z = omega_T2_rs_z * (current_rotation_difference - mTargetRotation);

                // Complete LHS contribution
                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_T2_r_x * omega_T2_s_x + second_order_term_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_T2_r_y * omega_T2_s_y + second_order_term_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_T2_r_z * omega_T2_s_z + second_order_term_z );
            }
        }

        // Then we consider the relation Slave-Slave ( SS )
        double divisor_s_1 = std::sqrt((1-omega_dot_T2_s*omega_dot_T2_s));
        double divisor_s_2 = std::pow((1-omega_dot_T2_s*omega_dot_T2_s),(3.0/2.0));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+2]);
            double omega_r_x_dot_T2 = inner_prod(omega_r_x,T2_s);
            double omega_r_y_dot_T2 = inner_prod(omega_r_y,T2_s);
            double omega_r_z_dot_T2 = inner_prod(omega_r_z,T2_s);
            double omega_T2_r_x = omega_r_x_dot_T2 / divisor_s_1;
            double omega_T2_r_y = omega_r_y_dot_T2 / divisor_s_1;
            double omega_T2_r_z = omega_r_z_dot_T2 / divisor_s_1;

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                Vector omega_s_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+0]);
                Vector omega_s_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+1]);
                Vector omega_s_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+2]);
                double omega_s_x_dot_T2 = inner_prod(omega_s_x,T2_s);
                double omega_s_y_dot_T2 = inner_prod(omega_s_y,T2_s);
                double omega_s_z_dot_T2 = inner_prod(omega_s_z,T2_s);
                double omega_T2_s_x = omega_s_x_dot_T2 / divisor_s_1;
                double omega_T2_s_y = omega_s_y_dot_T2 / divisor_s_1;
                double omega_T2_s_z = omega_s_z_dot_T2 / divisor_s_1;

                // First term of second derivative part
                Vector omega_rs_x = MathUtils<double>::CrossProduct(T3_s,t3rs_s[3*row_itr+0][3*collumn_itr+0]);
                Vector omega_rs_y = MathUtils<double>::CrossProduct(T3_s,t3rs_s[3*row_itr+1][3*collumn_itr+1]);
                Vector omega_rs_z = MathUtils<double>::CrossProduct(T3_s,t3rs_s[3*row_itr+2][3*collumn_itr+2]);

                double omega_T2_rs_x = inner_prod(omega_rs_x,T2_s) / divisor_s_1;
                double omega_T2_rs_y = inner_prod(omega_rs_y,T2_s) / divisor_s_1;
                double omega_T2_rs_z = inner_prod(omega_rs_z,T2_s) / divisor_s_1;

                // Second term of second derivative part
                omega_T2_rs_x += (omega_r_x_dot_T2 * omega_dot_T2_s * omega_s_x_dot_T2) / divisor_s_2;
                omega_T2_rs_y += (omega_r_y_dot_T2 * omega_dot_T2_s * omega_s_y_dot_T2) / divisor_s_2;
                omega_T2_rs_z += (omega_r_z_dot_T2 * omega_dot_T2_s * omega_s_z_dot_T2) / divisor_s_2;

                // Complete second derivative part
                double second_order_term_x = mSignFactor * omega_T2_rs_x * (current_rotation_difference - mTargetRotation);
                double second_order_term_y = mSignFactor * omega_T2_rs_y * (current_rotation_difference - mTargetRotation);
                double second_order_term_z = mSignFactor * omega_T2_rs_z * (current_rotation_difference - mTargetRotation);

                // Complete LHS contribution
                LHS( 3*row_id+0, 3*collumn_id+0 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_T2_r_x * omega_T2_s_x + second_order_term_x );
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_T2_r_y * omega_T2_s_y + second_order_term_y );
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mPenaltyFactor * mIntegrationWeight * mJ1 * ( omega_T2_r_z * omega_T2_s_z + second_order_term_z );
            }
        }

        // Then we consider the Master-Slave relation ( MS & SM )
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_T2_r_x = inner_prod(omega_r_x,T2_m) / divisor_m_1;
            double omega_T2_r_y = inner_prod(omega_r_y,T2_m) / divisor_m_1;
            double omega_T2_r_z = inner_prod(omega_r_z,T2_m) / divisor_m_1;

            for(int collumn_itr=0; collumn_itr<mNumberOfLocalEquationIdsOnSlave; collumn_itr++)
            {
                int collumn_id = mEquationIdsOfAffectedControlPointsOnSlave[collumn_itr];

                Vector omega_s_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+0]);
                Vector omega_s_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+1]);
                Vector omega_s_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*collumn_itr+2]);
                double omega_T2_s_x = inner_prod(omega_s_x,T2_s) / divisor_s_1;
                double omega_T2_s_y = inner_prod(omega_s_y,T2_s) / divisor_s_1;
                double omega_T2_s_z = inner_prod(omega_s_z,T2_s) / divisor_s_1;

                // MS
                LHS( 3*row_id+0, 3*collumn_id+0 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_r_x * omega_T2_s_x;
                LHS( 3*row_id+1, 3*collumn_id+1 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_r_y * omega_T2_s_y;
                LHS( 3*row_id+2, 3*collumn_id+2 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_r_z * omega_T2_s_z;

                // SM
                LHS( 3*collumn_id+0, 3*row_id+0 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_r_x * omega_T2_s_x;
                LHS( 3*collumn_id+1, 3*row_id+1 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_r_y * omega_T2_s_y;
                LHS( 3*collumn_id+2, 3*row_id+2 ) += mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * omega_T2_r_z * omega_T2_s_z;
            }
        }
    }

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
        double omega_dot_T2_m = inner_prod(omega_m, T2_m);
        double omega_dot_T2_s = inner_prod(omega_s, T2_s);
        double current_rotation_m = std::asin( omega_dot_T2_m );
        double current_rotation_s = std::asin( omega_dot_T2_s );
        double current_rotation_difference = current_rotation_m + mSignFactor*current_rotation_s;

        double target_rotation_in_deg = mTargetRotation*180/3.14;
        double current_rotation_difference_in_deg = current_rotation_difference*180/3.14;

        KRATOS_WATCH("--------------")
        // KRATOS_WATCH(w_m)
        // KRATOS_WATCH(w_s)
        // KRATOS_WATCH(omega_dot_T2_m)
        // KRATOS_WATCH(omega_dot_T2_s)
        KRATOS_WATCH(target_rotation_in_deg)
        KRATOS_WATCH(current_rotation_difference_in_deg)

        // Master contribution
        double divisor_m = std::sqrt((1-omega_dot_T2_m*omega_dot_T2_m));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnMaster; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnMaster[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*row_itr+2]);
            double omega_T2_r_x = inner_prod(omega_r_x,T2_m) / divisor_m;
            double omega_T2_r_y = inner_prod(omega_r_y,T2_m) / divisor_m;
            double omega_T2_r_z = inner_prod(omega_r_z,T2_m) / divisor_m;

            RHS[3*row_id+0] -= mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_T2_r_x;
            RHS[3*row_id+1] -= mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_T2_r_y;
            RHS[3*row_id+2] -= mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_T2_r_z;
        }

        // Slave contribution
        double divisor_s = std::sqrt((1-omega_dot_T2_s*omega_dot_T2_s));
        for(int row_itr=0; row_itr<mNumberOfLocalEquationIdsOnSlave; row_itr++)
        {
            int row_id = mEquationIdsOfAffectedControlPointsOnSlave[row_itr];

            Vector omega_r_x = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+0]);
            Vector omega_r_y = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+1]);
            Vector omega_r_z = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*row_itr+2]);
            double omega_T2_r_x = inner_prod(omega_r_x,T2_s) / divisor_s;
            double omega_T2_r_y = inner_prod(omega_r_y,T2_s) / divisor_s;
            double omega_T2_r_z = inner_prod(omega_r_z,T2_s) / divisor_s;

            RHS[3*row_id+0] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_T2_r_x;
            RHS[3*row_id+1] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_T2_r_y;
            RHS[3*row_id+2] -= mSignFactor * mPenaltyFactor * mIntegrationWeight * mJ1 * (current_rotation_difference - mTargetRotation) * omega_T2_r_z;
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
