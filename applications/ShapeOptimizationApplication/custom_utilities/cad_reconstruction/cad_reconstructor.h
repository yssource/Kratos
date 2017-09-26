// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_RECONSTRUCTOR_H
#define CAD_RECONSTRUCTOR_H

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
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"

// ==============================================================================

namespace Kratos
{
class CADReconstructor
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CADReconstructor
    KRATOS_CLASS_POINTER_DEFINITION(CADReconstructor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADReconstructor()
    {      
    }

    /// Destructor.
    virtual ~CADReconstructor()
    {
    }

    // --------------------------------------------------------------------------
    // void CreateReconstructionDataBase(ModelPart& fe_model_part, boost::python::dict cad_geometry, boost::python::dict cad_integration_data)
    // {
    //   mReconstructionDataBase = ReconstructionDataBase(fe_model_part, cad_geometry, cad_integration_data);
    // }

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "CADReconstructor";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "CADReconstructor";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

  private:

	// ==============================================================================
    // Solver and strategies
    // ==============================================================================
	// CompressedLinearSolverType& m_linear_solver;

    /// Assignment operator.
    //      CADReconstructor& operator=(CADReconstructor const& rOther);

    /// Copy constructor.
    //      CADReconstructor(CADReconstructor const& rOther);

}; // Class CADReconstructor
} // namespace Kratos.

#endif // CAD_RECONSTRUCTOR_H
