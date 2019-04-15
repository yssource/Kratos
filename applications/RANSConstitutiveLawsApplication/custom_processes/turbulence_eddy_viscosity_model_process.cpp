#include "turbulence_eddy_viscosity_model_process.h"

#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
/* Public functions *******************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::TurbulenceEddyViscosityModelProcess(
    ModelPart& rModelPart, Process& rRansYPlusProcess)
    : mrModelPart(rModelPart), mrRansYPlusProcess(rRansYPlusProcess)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Execute()
{
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
std::string TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Info() const
{
    return "TurbulenceEddyViscosityModelProcess";
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "TurbulenceEddyViscosityModelProcess";
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PrintData(
    std::ostream& rOStream) const
{
}

/* Protected functions ****************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::CalculateYplus(unsigned int Step)
{
    mrRansYPlusProcess.Execute();
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::UpdateFluidViscosity()
{
    KRATOS_TRY

    NodesArrayType& nodes = mrModelPart.Nodes();

// Modifying viscosity of the nodes with the calculated turbulent viscosity
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
    {
        auto it_node = nodes.begin() + i;
        const double kinematic_viscosity =
            it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double turbulent_viscosity =
            it_node->FastGetSolutionStepValue(TURBULENT_VISCOSITY);

        double& effective_viscosity = it_node->FastGetSolutionStepValue(VISCOSITY);
        effective_viscosity = kinematic_viscosity + turbulent_viscosity;
    }

    KRATOS_CATCH("");
}

/* Protected functions ****************************************************/

/* External functions *****************************************************/

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
typedef LinearSolver<SparseSpaceType, DenseSpaceType> LinearSolverType;

template class TurbulenceEddyViscosityModelProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType>;
template class TurbulenceEddyViscosityModelProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType>;
} // namespace Kratos
