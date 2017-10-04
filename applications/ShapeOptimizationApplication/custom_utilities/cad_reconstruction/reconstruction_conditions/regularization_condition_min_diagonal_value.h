// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef REGULARIZATION_CONDITION_MIN_DIAGONAL_VALUE_H
#define REGULARIZATION_CONDITION_MIN_DIAGONAL_VALUE_H

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
#include "regularization_condition_base.h"

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

class MinimalDiagonalValueCondition : public RegularizationCondition
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Pointer definition of MinimalDiagonalValueCondition
    KRATOS_CLASS_POINTER_DEFINITION(MinimalDiagonalValueCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MinimalDiagonalValueCondition( double min_value )
    : mMinValue( min_value )
    {
    }

    /// Destructor.
    virtual ~MinimalDiagonalValueCondition()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void Initialize()
    {
    }  
    
    // --------------------------------------------------------------------------
    void ComputeAndAddLHSContribution( CompressedMatrix& LHS )
    {
        std::cout << "> Starting regularization of LHS by enforcing minimal diagonal value..." << std::endl;           
    
        double zero_threshold = 1e-10;

        for(int i_itr=0; i_itr<LHS.size1(); i_itr++)
            if( std::abs(LHS(i_itr,i_itr)) < zero_threshold )
                LHS(i_itr,i_itr) = mMinValue;

        std::cout << "> All values on main diagonal of LHS < " << zero_threshold << " are manually set to " << mMinValue << "." << std::endl;
        std::cout << "> Finished regularization of LHS by enforcing minimal diagonal value." << std::endl;           
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS )
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
        return "MinimalDiagonalValueCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "MinimalDiagonalValueCondition";
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    double mMinValue;

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
    //      MinimalDiagonalValueCondition& operator=(MinimalDiagonalValueCondition const& rOther);

    /// Copy constructor.
    //      MinimalDiagonalValueCondition(MinimalDiagonalValueCondition const& rOther);

    ///@}

}; // Class MinimalDiagonalValueCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // REGULARIZATION_CONDITION_MIN_DIAGONAL_VALUE_H
