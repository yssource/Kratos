// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef REGULARIZATION_CONDITION_MIN_CP_DISPLACEMENT_H
#define REGULARIZATION_CONDITION_MIN_CP_DISPLACEMENT_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

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

class MinimalControlPointDisplacementCondition : public RegularizationCondition
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Pointer definition of MinimalControlPointDisplacementCondition
    KRATOS_CLASS_POINTER_DEFINITION(MinimalControlPointDisplacementCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MinimalControlPointDisplacementCondition( ReconstructionDataBase& data_base, double beta_value, std::string solution_strategy )
    : mrReconstructionDataBase( data_base ),
      mBetaValue( beta_value ),
      mSolutionStrategy( solution_strategy )
    {
    }

    /// Destructor.
    virtual ~MinimalControlPointDisplacementCondition()
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
        std::cout << "> Starting to compute LHS contribution of beta-regularization (minimal control point displacement)..." << std::endl;  

        for(size_t i_itr=0; i_itr<LHS.size1(); i_itr++)
            LHS(i_itr,i_itr) += mBetaValue;
            
        std::cout << "> Finished computing LHS contribution of beta-regularization (minimal control point displacement)." << std::endl;           
    }

    // --------------------------------------------------------------------------
    void ComputeAndAddRHSContribution( Vector& RHS )
    {
        std::cout << "> Starting to compute RHS contribution of beta-regularization (minimal control point displacement)..." << std::endl;  

        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
            {
                if(control_point_i.IsRelevantForReconstruction())
                {
                    unsigned int cp_equation_id = control_point_i.GetEquationId();
                    RHS[3*cp_equation_id+0] -= mBetaValue*control_point_i.GetdX();
                    RHS[3*cp_equation_id+1] -= mBetaValue*control_point_i.GetdY();
                    RHS[3*cp_equation_id+2] -= mBetaValue*control_point_i.GetdZ();     
                }
            }
        }
        
        std::cout << "> Finished computing RHS contribution of beta-regularization (minimal control point displacement)." << std::endl;                   
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
        return "MinimalControlPointDisplacementCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "MinimalControlPointDisplacementCondition";
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
    ReconstructionDataBase& mrReconstructionDataBase;
    double mBetaValue;
    std::string mSolutionStrategy;

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
    //      MinimalControlPointDisplacementCondition& operator=(MinimalControlPointDisplacementCondition const& rOther);

    /// Copy constructor.
    //      MinimalControlPointDisplacementCondition(MinimalControlPointDisplacementCondition const& rOther);

    ///@}

}; // Class MinimalControlPointDisplacementCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // REGULARIZATION_CONDITION_MIN_CP_DISPLACEMENT_H
