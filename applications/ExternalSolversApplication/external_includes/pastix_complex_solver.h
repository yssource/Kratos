//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:         BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $ExternalSolversApplication   $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:             April 2017 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_PASTIX_COMPLEX_SOLVER)
#define KRATOS_PASTIX_COMPLEX_SOLVER

// System includes
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <complex>
#include <vector>
#include <algorithm>

// External includes
extern "C" {
#include <pastix.h>
}

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class PastixComplexSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PastixComplexSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    PastixComplexSolver(Parameters& r_settings)
    {
        mp_pastix_data = NULL;

        Parameters default_settings(R"(
            {
                "solver_type" : "pastix",
                "echo_level" : 0
            })");

        r_settings.ValidateAndAssignDefaults(default_settings);

        m_echo_level = r_settings["echo_level"].GetInt();
        if (m_echo_level > 4 || m_echo_level < 0)
            KRATOS_ERROR << "invalid echo_level: " << m_echo_level << std::endl;
    }

    PastixComplexSolver(const PastixComplexSolver& Other) = delete;

    ~PastixComplexSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    PastixComplexSolver& operator=(const PastixComplexSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mp_pastix_data != NULL)
            this->Clear();

        m_nrows = rA.size1();

        if (m_rowptr.size() != rA.index1_data().size())
            m_rowptr.resize(rA.index1_data().size());

        if (m_col.size() != rA.index2_data().size())
            m_col.resize(rA.index2_data().size());

        if (m_perm.size() != static_cast<unsigned int>(m_nrows))
            m_perm.resize(m_nrows);

        if (m_invp.size() != static_cast<unsigned int>(m_nrows))
            m_invp.resize(m_nrows);

        // initialize iparm and dparm values
#ifdef VERSION_PASTIX_6
        m_iparm[IPARM_MODIFY_PARAMETER] = 0;
        m_iparm[IPARM_START_TASK      ] = PastixTaskInit;
        m_iparm[IPARM_END_TASK        ] = PastixTaskInit;
#else
        m_iparm[IPARM_MODIFY_PARAMETER] = API_NO;
        m_iparm[IPARM_START_TASK      ] = API_TASK_INIT;
        m_iparm[IPARM_END_TASK        ] = API_TASK_INIT;
#endif

        z_pastix(&mp_pastix_data, 0, m_nrows, NULL, NULL, NULL, NULL, NULL, NULL, 1, m_iparm, m_dparm);
    }

    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mp_pastix_data != NULL)
            KRATOS_ERROR << "pastix_data != NULL upon entering InitializeSolutionStep." << std::endl;

        // 1-based indexing for pastix
        std::copy(std::begin(rA.index1_data()),std::end(rA.index1_data()),std::begin(m_rowptr));
        for (auto &i : m_rowptr)
            ++i;
        std::copy(std::begin(rA.index2_data()),std::end(rA.index2_data()),std::begin(m_col));
        for (auto &i : m_col)
            ++i;

        // Factorize system matrix
#ifdef VERSION_PASTIX_6
        m_iparm[IPARM_VERBOSE        ] = static_cast<pastix_verbose_e>(m_echo_level);
        m_iparm[IPARM_FACTORIZATION  ] = PastixFactLU; // LU factorization
#else
        m_iparm[IPARM_VERBOSE        ] = static_cast<API_VERBOSE>(m_echo_level);
        m_iparm[IPARM_RHS_MAKING     ] = API_RHS_B; // user-provided rhs
        m_iparm[IPARM_SYM            ] = API_SYM_NO; // non-symmetric
        m_iparm[IPARM_FACTORIZATION  ] = API_FACT_LU; // LU factorization
        m_iparm[IPARM_TRANSPOSE_SOLVE] = API_YES; // solve transpose for csr matrix format
#endif
#ifdef _OPENMP
        m_iparm[IPARM_THREAD_NBR     ] = omp_get_max_threads();
#endif

#ifdef VERSION_PASTIX_6
        m_iparm[IPARM_START_TASK] = PastixTaskOrdering;
        m_iparm[IPARM_END_TASK  ] = PastixTaskNumfact;
#else
        m_iparm[IPARM_START_TASK] = API_TASK_ORDERING;
        m_iparm[IPARM_END_TASK  ] = API_TASK_NUMFACT;
#endif

        z_pastix(&mp_pastix_data, 0, m_nrows,
            static_cast<pastix_int_t*>(&m_rowptr[0]),
            static_cast<pastix_int_t*>(&m_col[0]),
            static_cast<std::complex<double>*>(&rA.value_data()[0]),
            static_cast<pastix_int_t*>(&m_perm[0]),
            static_cast<pastix_int_t*>(&m_invp[0]),
            NULL,
            1, m_iparm, m_dparm);
    }

    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mp_pastix_data == NULL)
            KRATOS_ERROR << "pastix_data == NULL upon entering Solve." << std::endl;

        noalias(rX) = rB;
#ifdef VERSION_PASTIX_6
        m_iparm[IPARM_START_TASK] = PastixTaskSolve;
        m_iparm[IPARM_END_TASK  ] = PastixTaskSolve;
#else
        m_iparm[IPARM_START_TASK] = API_TASK_SOLVE;
        m_iparm[IPARM_END_TASK  ] = API_TASK_SOLVE;
#endif

        z_pastix(&mp_pastix_data, 0, m_nrows,
            static_cast<pastix_int_t*>(&m_rowptr[0]),
            static_cast<pastix_int_t*>(&m_col[0]),
            static_cast<std::complex<double>*>(&rA.value_data()[0]),
            static_cast<pastix_int_t*>(&m_perm[0]),
            static_cast<pastix_int_t*>(&m_invp[0]),
            static_cast<std::complex<double>*>(&rX[0]),
            1, m_iparm, m_dparm);
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);
        FinalizeSolutionStep(rA, rX, rB);

        return true;
    }

    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        this->Clear();
    }

    void Clear() override
    {
        if (mp_pastix_data != NULL)
        {
#ifdef VERSION_PASTIX_6
            m_iparm[IPARM_START_TASK] = PastixTaskClean;
            m_iparm[IPARM_END_TASK  ] = PastixTaskClean;
#else
            m_iparm[IPARM_START_TASK] = API_TASK_CLEAN;
            m_iparm[IPARM_END_TASK  ] = API_TASK_CLEAN;
#endif

            z_pastix(&mp_pastix_data, 0, m_nrows,
                static_cast<pastix_int_t*>(&m_rowptr[0]),
                static_cast<pastix_int_t*>(&m_col[0]),
                NULL,
                static_cast<pastix_int_t*>(&m_perm[0]),
                static_cast<pastix_int_t*>(&m_invp[0]),
                NULL,
                1, m_iparm, m_dparm);

            mp_pastix_data = NULL;
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Pastix direct solver finished.";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    pastix_data_t* mp_pastix_data;
    pastix_int_t m_nrows;
    std::vector<pastix_int_t> m_rowptr;
    std::vector<pastix_int_t> m_col;
    std::vector<pastix_int_t> m_perm;
    std::vector<pastix_int_t> m_invp;
    pastix_int_t m_iparm[IPARM_SIZE];
    double m_dparm[DPARM_SIZE];
    int m_echo_level;

    ///@}
}; // Class PastixComplexSolver

///@}

///@name Input and output
///@{

/**
 * input stream function
 */
template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator>>(std::istream& rIStream,
                                PastixComplexSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const PastixComplexSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_PASTIX_COMPLEX_SOLVER defined
