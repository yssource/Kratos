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
    void ComputeAndAddLHSContribution( SparseMatrix& LHS ) override
    {
                // return mrAffectedPatch.GetEquationIdsOfAffectedControlPoints( mParmeterSpans, mParmeterValues );    
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS ) override
    {
    }

    // --------------------------------------------------------------------------
    void FlagControlPointsRelevantForReconstruction()
    {
        mrAffectedPatch.FlagAffectedControlPointsForReconstruction(  mParmeterSpans, mParmeterValues );
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
