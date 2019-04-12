//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED)
#define K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes
#include "custom_processes/scalar_co_solving_process.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// The base class for all KEpsilonCoSolvingProcess in Kratos.
/** The KEpsilonCoSolvingProcess is the base class for all KEpsilonCoSolvingProcess and defines a simple interface for them.
    Execute method is used to execute the KEpsilonCoSolvingProcess algorithms. While the parameters of this method
  can be very different from one KEpsilonCoSolvingProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all KEpsilonCoSolvingProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other KEpsilonCoSolvingProcess or the base KEpsilonCoSolvingProcess class.
*/
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class KEpsilonCoSolvingProcess
    : public ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    typedef ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// Pointer definition of KEpsilonCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(KEpsilonCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    KEpsilonCoSolvingProcess(ModelPart& rModelPart, Parameters& rParameters, Process& rYPlusModelProcess)
        : BaseType(rModelPart, rParameters, TURBULENT_VISCOSITY), mrYPlusModelProcess(rYPlusModelProcess)
    {
    }

    /// Destructor.
    ~KEpsilonCoSolvingProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
    std::string Info() const override
    {
        return "KEpsilonCoSolvingProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KEpsilonCoSolvingProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Process& mrYPlusModelProcess;

    ///@}
    ///@name Operations
    ///@{

    void UpdateBeforeSolveSolutionStep() override
    {
        mrYPlusModelProcess.Execute();
    }

    void UpdateAfterSolveSolutionStep() override
    {
        EvmKepsilonModelUtilities().CalculateTurbulentViscosityForModelPart(this->mrModelPart);
    }

    void UpdateConvergenceVariable() override
    {
        RansCalculationUtilities().UpdateEffectiveViscosityForModelPart(this->mrModelPart);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(
        KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    // KEpsilonCoSolvingProcess(KEpsilonCoSolvingProcess const& rOther);

    ///@}

}; // Class KEpsilonCoSolvingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::istream& operator>>(std::istream& rIStream,
                                KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED defined
