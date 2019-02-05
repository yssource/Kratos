#include "turbulence_evm_k_epsilon_process.h"

namespace Kratos
{
/* Public functions *******************************************************/
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::TurbulenceEvmKEpsilonProcess(
    ModelPart& rModelPart, Parameters& rParameters, TLinearSolver& rLinearSolver)
    : Process(), mrModelPart(rModelPart), mrParameters(rParameters), mrLinearSolver(rLinearSolver)
{

}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::~TurbulenceEvmKEpsilonProcess()
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::Execute()
{
    this->AssignBoussinesqForce();
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitialize()
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitializeSolutionStep()
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
std::string TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::Info() const
{
    return "TurbulenceEvmKEpsilonProcess";
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "TurbulenceEvmKEpsilonProcess";
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::PrintData(
    std::ostream& rOStream) const
{
}

/* Protected functions ****************************************************/

/* External functions *****************************************************/

/// output stream function
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos
