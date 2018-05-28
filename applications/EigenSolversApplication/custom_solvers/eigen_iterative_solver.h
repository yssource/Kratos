/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_EIGEN_ITERATIVE_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_ITERATIVE_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_utilities/ublas_wrapper.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "eigen_direct_solver.h" // for the typenames

namespace Kratos
{

template <typename scalar_t>
struct ConjugateGradient : public SolverType<scalar_t>
{
    using TSolver = Eigen::ConjugateGradient<typename SolverType<scalar_t>::TSparseMatrix,
                                            Eigen::Lower|Eigen::Upper >;

    static constexpr auto Name = "ConjugateGradient";
};

template <typename scalar_t>
struct BiCGSTAB : public SolverType<scalar_t>
{
    using TSolver = Eigen::BiCGSTAB<typename SolverType<scalar_t>::TSparseMatrix >;

    static constexpr auto Name = "BiCGSTAB";
};

template <
    class TSolverType,
    class TSparseSpaceType = typename TSolverType::TGlobalSpace,
    class TDenseSpaceType = typename TSolverType::TLocalSpace>
class EigenIterativeSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType>
{
    Parameters mParam;

    typename TSolverType::TSolver m_solver;

    EigenIterativeSolver &operator=(const EigenIterativeSolver &Other);

    EigenIterativeSolver(const EigenIterativeSolver &Other);

  public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenIterativeSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigenIterativeSolver(Parameters param) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "eigen_cg; eigen_bicgstab",
            "max_iteration": 1000,
            "tolerance": 1e-8,
            "echo_level": 1
        })");

        // Eigen uses n_dof as max iterations

        mParam.ValidateAndAssignDefaults(default_params);
    }

    ~EigenIterativeSolver() override {}

    /**
     * This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * For example if we are implementing a Iterative solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {

    }

    /**
     * This function actually performs the solution work, eventually taking advantage of what was done before in the
     * Initialize and InitializeSolutionStep functions.
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     */
    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        const int echo_level = mParam["echo_level"].GetInt();

        // TODO: wrapper connot be deleted before solution is completed, because it stores the indices as a member.
        // for the iterative solver, they are still needed.
        UblasWrapper<DataType> ublas_wrapper(rA);

        const auto& a = ublas_wrapper.matrix();

        m_solver.setMaxIterations(mParam["max_iteration"].GetInt());
        m_solver.setTolerance(mParam["tolerance"].GetDouble());

        m_solver.compute(a);

        KRATOS_ERROR_IF(m_solver.info() != Eigen::Success) << "InitializeSolutionStep failed!" << std::endl;

        if (echo_level > 0) {
            std::cout << "Initialized iterative solver."  << std::endl;
        }

        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1> > x(rX.data().begin(), rX.size());
        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1> > b(rB.data().begin(), rB.size());

        x = m_solver.solveWithGuess(b, x);

        if (echo_level > 0) {
            std::cout << "Completed iterative solver."  << std::endl;
            std::cout << "#iterations:     " << m_solver.iterations() << std::endl;
            std::cout << "estimated error: " << m_solver.error()      << std::endl;
        }

        KRATOS_ERROR_IF(m_solver.info() != Eigen::Success) << "Solving failed!" << std::endl;
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    bool Solve(SparseMatrixType &rA, VectorType &rX, VectorType &rB) override
    {
        InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);

        return true;
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution matrix
     * @param rB Right hand side matrix
     * @return true if solution found, otherwise false
     */
    bool Solve(SparseMatrixType &rA, DenseMatrixType &rX, DenseMatrixType &rB) override
    {
        return false;
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "EigenIterativeSolver<" << TSolverType::Name << "> finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }
}; // class EigenIterativeSolver

/**
 * input stream function
 */
template<
    class TSolverType,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType>
inline std::istream &operator>>(
    std::istream &rIStream,
    EigenIterativeSolver<TSolverType, TSparseSpaceType, TDenseSpaceType> &rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<
    class TSolverType,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType
    >
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EigenIterativeSolver<TSolverType, TSparseSpaceType, TDenseSpaceType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_ITERATIVE_SOLVER_H_INCLUDED)
