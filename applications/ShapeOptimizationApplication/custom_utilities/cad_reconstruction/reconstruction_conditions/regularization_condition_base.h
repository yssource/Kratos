// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef REGULARIZATION_CONDITION_H
#define REGULARIZATION_CONDITION_H

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
#include "includes/define.h"

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

class RegularizationCondition
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Pointer definition of RegularizationCondition
    KRATOS_CLASS_POINTER_DEFINITION(RegularizationCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RegularizationCondition()
    {
    }

    /// Destructor.
    virtual ~RegularizationCondition()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    virtual void FlagControlPointsRelevantForReconstruction(){};
    
    // --------------------------------------------------------------------------    
    virtual void Initialize() = 0;    
    
    // --------------------------------------------------------------------------
    virtual void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) = 0;

    // --------------------------------------------------------------------------
    virtual void ComputeAndAddRHSContribution( Vector& RHS ) = 0;

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
        return "RegularizationCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "RegularizationCondition";
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
    //      RegularizationCondition& operator=(RegularizationCondition const& rOther);

    /// Copy constructor.
    //      RegularizationCondition(RegularizationCondition const& rOther);

    ///@}

}; // Class RegularizationCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // REGULARIZATION_CONDITION_H
