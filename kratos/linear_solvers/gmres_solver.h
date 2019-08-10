//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Denis Demidov
//    Adapted by:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_AMGCL_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_AMGCL_GMRES_SOLVER_H_INCLUDED

// System includes

// External includes
/* AMGCL */
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>

// #ifdef AMGCL_GPGPU // TODO: Include GPU support
// #  include <amgcl/backend/vexcl.hpp>
// #  include <amgcl/backend/vexcl_static_matrix.hpp>
// #endif

// Project includes
#include "includes/ublas_interface.h"
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"

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

/**
 * @class AMGCLGMRESSolver
 * @ingroup KratosCore
 * @brief This is an adaptation of the GMRES solver included in the AMGCL library
 * @details This implementation can be extended the other solvers available in AMGCL
 * @tparam TSparseSpaceType The sparse matrix type
 * @tparam TDenseSpaceType The dense matrix type
 * @tparam TPreconditionerType The preconditioner type
 * @tparam TReordererType The reorder type
 * @author Denis Demidov
 * @author Vicente Mataix Ferrandiz
 * @todo Clean up old versions of GMRES in Kratos code
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AMGCLGMRESSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AMGCLGMRESSolver
    KRATOS_CLASS_POINTER_DEFINITION(  AMGCLGMRESSolver );

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AMGCLGMRESSolver() {}

    /**
     * @brief Constructor with tolerance
     * @param MaxTolerance The tolerance of convergence
     */
    AMGCLGMRESSolver(const double MaxTolerance) : BaseType(MaxTolerance) {}

    /**
     * @brief Constructor with tolerance and max number of iterations
     * @param MaxTolerance The tolerance of convergence
     * @param MaxIterationsNumber The max number of iterations
     */
    AMGCLGMRESSolver(
        const double MaxTolerance,
        const std::size_t MaxIterationsNumber
        ) : BaseType(MaxTolerance, MaxIterationsNumber)
    {

    }

    /**
     * @brief Constructor with tolerance, max number of iterations and the preconditioner pointer
     * @param MaxTolerance The tolerance of convergence
     * @param MaxIterationsNumber The max number of iterations
     * @param pPreconditioner The preconditioner considered
     */
    AMGCLGMRESSolver(
        const double MaxTolerance,
        const std::size_t MaxIterationsNumber,
        typename TPreconditionerType::Pointer pPreconditioner
        ) : BaseType(MaxTolerance, MaxIterationsNumber, pPreconditioner)
    {

    }

    /**
     * @brief Constructor with parameters
     * @param SettingsParameters The parameters considered
     */
     AMGCLGMRESSolver(Parameters SettingsParameters)
        : BaseType(SettingsParameters),
          mSettingsParameters(SettingsParameters)
    {
        Parameters default_parameters = GetDefaultParameters();
        mSettingsParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }

    /**
     * @brief Constructor with parameters and the preconditioner pointer
     * @param SettingsParameters The parameters considered
     * @param pPreconditioner The preconditioner considered
     */
     AMGCLGMRESSolver(
         Parameters SettingsParameters,
         typename TPreconditionerType::Pointer pPreconditioner
        ) : BaseType(SettingsParameters, pPreconditioner),
            mSettingsParameters(SettingsParameters)
    {
        Parameters default_parameters = GetDefaultParameters();
        mSettingsParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }

    /// Copy constructor.
    AMGCLGMRESSolver(const AMGCLGMRESSolver& Other)
        : BaseType(Other)
    {

    }

    /// Destructor.
    ~AMGCLGMRESSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AMGCLGMRESSolver& operator=(const AMGCLGMRESSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Normal solve method.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        bool is_solved= false;


        return is_solved;
    }

    /** @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    bool Solve(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB
        ) override
    {
        //DOES NOTHING AT THE MOMENT
        bool is_solved = true;


        return is_solved;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "GMRES iterative solver [at the moment unfortunately without] with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << "GMRES iterative solver [at the moment unfortunately without] with ";
        BaseType::GetPreconditioner()->PrintInfo(OStream);
    }

    /// Print object's data.
    void  PrintData(std::ostream& OStream) const override
    {
        OStream << "GMRES configuration parameters: " << mSettingsParameters << std::endl;
        BaseType::PrintData(OStream);
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

    Parameters mSettingsParameters;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {

        })" );

        return default_parameters;
    }

    ///@}

}; // Class AMGCLGMRESSolver
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& rIStream,
                                  AMGCLGMRESSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const AMGCLGMRESSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
//   ///@}
//
//
}  // namespace Kratos.

#endif // KRATOS_AMGCL_GMRES_SOLVER_H_INCLUDED  defined


