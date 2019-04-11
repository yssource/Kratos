// System includes

#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// Application includes
// #include "custom_processes/turbulence_eddy_viscosity_model_process.h"
// #include "custom_processes/turbulence_evm_k_epsilon_process.h"
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_process.h"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    // typedef TurbulenceEvmKEpsilonProcess<2, SparseSpaceType, LocalSpaceType, LinearSolverType> TurbulenceEvmKEpsilon2DProcess;
    // py::class_<TurbulenceEvmKEpsilon2DProcess, TurbulenceEvmKEpsilon2DProcess::Pointer, Process>(
    //     m, "TurbulenceEvmKEpsilon2DProcess")
    //     .def(py::init<ModelPart&, Parameters&, LinearSolverType::Pointer,
    //                   LinearSolverType::Pointer, LinearSolverType::Pointer>());

    // typedef TurbulenceEvmKEpsilonProcess<3, SparseSpaceType, LocalSpaceType, LinearSolverType> TurbulenceEvmKEpsilon3DProcess;
    // py::class_<TurbulenceEvmKEpsilon3DProcess, TurbulenceEvmKEpsilon3DProcess::Pointer, Process>(
    //     m, "TurbulenceEvmKEpsilon3DProcess")
    //     .def(py::init<ModelPart&, Parameters&, LinearSolverType::Pointer,
    //                   LinearSolverType::Pointer, LinearSolverType::Pointer>());

    py::class_<RansLogarithmicYPlusModelProcess, RansLogarithmicYPlusModelProcess::Pointer, Process>(
        m, "RansLogarithmicYPlusModelProcess")
        .def(py::init<ModelPart&, Parameters&>());
}

} // namespace Python.
} // Namespace Kratos
