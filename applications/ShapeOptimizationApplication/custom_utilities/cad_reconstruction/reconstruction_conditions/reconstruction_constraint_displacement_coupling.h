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
    DisplacementCouplingConstraint( double penalty_factor )
    : mPenaltyFactor( penalty_factor )
    {
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
    }

    // --------------------------------------------------------------------------
    void Initialize() override
    {
    }

    // --------------------------------------------------------------------------    
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) override
    {
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
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
    double mPenaltyFactor = 0;
    
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
