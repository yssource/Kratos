//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein
//

#if !defined(KRATOS_AMESOS2_SOLVER_H_INCLUDED )
#define  KRATOS_AMESOS2_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "trilinos_space.h"
#include "Epetra_FEVector.h"

//amesos2 solver includes
#include "Teuchos_RCP.hpp"
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#include "Epetra_LinearProblem.h"



namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Amesos2Solver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of AmesosSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(Amesos2Solver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> MPISparseSpaceType;

    Amesos2Solver(Parameters settings)
    {
        Parameters default_settings( R"(
        {
        "solver_type": "Amesos2Solver",
        "amesos2_solver_type" : "Klu2",
        "trilinos_amesos2_parameter_list": {
            }
        }  )" );

        settings.ValidateAndAssignDefaults(default_settings);

        //assign the amesos2 parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mParameterList
        mParameterList = Teuchos::ParameterList();
        for(auto it = settings["trilinos_amesos2_parameter_list"].begin(); it != settings["trilinos_amesos2_parameter_list"].end(); it++)
        {
            if(it->IsString()) mParameterList.set(it.name(), it->GetString());
            else if(it->IsInt()) mParameterList.set(it.name(), it->GetInt());
            else if(it->IsBool()) mParameterList.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mParameterList.set(it.name(), it->GetDouble());
        }

        mSolverName = settings["amesos2_solver_type"].GetString();

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos2 solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /**
     * Default constructor
     */
    Amesos2Solver(const std::string& SolverName, Teuchos::ParameterList& rParameterList)
    {
        mParameterList = rParameterList;
        mSolverName = SolverName;


        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos2 solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;

    }

    /**
     * Destructor
     */
    virtual ~Amesos2Solver() {}
    
   
    static bool HasSolver(const std::string& Amesos2SolverName)
    {
        return Amesos2::query(Amesos2SolverName);
    }
    

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        KRATOS_TRY
        rA.Comm().Barrier();
         
	    Amesos2::Solver<Epetra_FECrsMatrix, Epetra_FEVector> *p_amesos2_solver = Amesos2::create<Epetra_FECrsMatrix, Epetra_FEVector> (mSolverName, &rA, &rX, &rB);

        //p_amesos2_solver->SetParameters( mParameterList );
        //p_amesos2_solver->SymbolicFactorization();
        //p_amesos2_solver->NumericFactorization();
        //p_amesos2_solver->Solve();

        //delete p_amesos2_solver;

        rA.Comm().Barrier();

        return true;
        KRATOS_CATCH("");
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Amesos2 solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    Teuchos::ParameterList mParameterList;
    std::string mSolverName;

    /**
     * Assignment operator.
     */
    Amesos2Solver& operator=(const Amesos2Solver& Other);

    /**
     * Copy constructor.
     */
    Amesos2Solver(const Amesos2Solver& Other);

}; // Class Amesos2Solver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, Amesos2Solver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Amesos2Solver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_AMESOS2_SOLVER_H_INCLUDED  defined

