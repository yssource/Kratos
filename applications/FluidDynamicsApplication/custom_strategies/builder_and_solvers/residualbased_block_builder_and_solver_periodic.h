//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#ifndef KRATOS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H
#define KRATOS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H

/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Variant of ResidualBasedBlockBuilderAndSolver for problems with periodic boundary conditions.
 * @see PeriodicCondition
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolverPeriodic
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverPeriodic);


    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedBlockBuilderAndSolverPeriodic(typename TLinearSolver::Pointer pNewLinearSystemSolver,
                                               const Variable<int>& PeriodicVariable ):
        ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver),
        mPeriodicIdVar(PeriodicVariable)
    {}

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverPeriodic() override
    {
    }


    /*@} */
    /**@name Operators */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */


    void SetUpSystem(ModelPart &r_model_part) override
    {
        // Assign an Equation Id to all non-duplicate nodes
        unsigned int EqId = 0;

        for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin(); itDof != BaseType::mDofSet.end(); ++itDof)
        {
            if ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id()) ) {
                itDof->SetEquationId(EqId++);
            }
            else {

            }
        }

        // Copy Equation Id to duplicate nodes.
        for (ModelPart::ConditionIterator itCond = r_model_part.ConditionsBegin(); itCond != r_model_part.ConditionsEnd(); itCond++)
        {
            // PeriodicCondition always have exactly 2 nodes
            ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
            if ( (itCond->Is(PERIODIC)) && (rGeom.PointsNumber() == 2) )
            {
                int Node0 = rGeom[0].Id();
                int Node0Pair = rGeom[0].FastGetSolutionStepValue(mPeriodicIdVar);

                int Node1 = rGeom[1].Id();
                int Node1Pair = rGeom[1].FastGetSolutionStepValue(mPeriodicIdVar);

                // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                if ( ( Node0 == Node1Pair ) && ( Node1 == Node0Pair ) )
                {
                    KRATOS_ERROR_IF_NOT(itCond->GetProperties().Has(PERIODIC_VARIABLES)) <<
                        "Condition " << itCond->Id() << " is marked as periodic but does not have PERIODIC_VARIABLES defined in its properties." << std::endl;

                    const auto& periodic_variable_list = itCond->GetProperties().GetValue(PERIODIC_VARIABLES);

                    // The condition checks if Node0 is the one with lower Id of the pair (that is: the one that does not have an EquationId yet)
                    Node<3>& r_origin_node = ( Node0 < Node0Pair ) ? rGeom[1] : rGeom[0];
                    Node<3>& r_destination_node = ( Node0 < Node0Pair ) ? rGeom[0] : rGeom[1];

                    for (auto var = periodic_variable_list.DoubleVariablesBegin(); var != periodic_variable_list.DoubleVariablesEnd(); ++var) {
                        CopyEquationId(*var, r_origin_node, r_destination_node);
                    }

                    for (auto var = periodic_variable_list.VariableComponentsBegin(); var != periodic_variable_list.VariableComponentsEnd(); ++var) {
                        CopyEquationId(*var, r_origin_node, r_destination_node);
                    }

                    // Assign an Equation Id for any remaining DOF in the destination node
                    for ( Node<3>::DofsContainerType::iterator i_dof = r_destination_node.GetDofs().begin(); 
                            i_dof != r_destination_node.GetDofs().end(); ++i_dof) {
                        if (i_dof->EquationId() != r_origin_node.pGetDof(i_dof->GetVariable())->EquationId()) {
                            i_dof->SetEquationId(EqId++);
                        }
                    }
                }
            }
        }

        BaseType::mEquationSystemSize = EqId;
    }

    void ApplyDirichletConditions(typename TSchemeType::Pointer pScheme,
                                          ModelPart &r_model_part,
                                          TSystemMatrixType &A,
                                          TSystemVectorType &Dx,
                                          TSystemVectorType &b) override
    {
        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin(); itDof != BaseType::mDofSet.end(); ++itDof)
        {
            if (itDof->IsFixed())
            {
                std::size_t RowId = itDof->EquationId();
                std::size_t RowBegin = Arow_indices[RowId];
                std::size_t RowEnd = Arow_indices[RowId+1];

                for (std::size_t k = RowBegin; k != RowEnd; k++)
                {
                    if ( Acol_indices[k] == RowId )
                        Avalues[k] = 1.0;
                    else
                        Avalues[k] = 0.0;
                }

                b[RowId] = 0.0;
            }
        }
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    const Variable<int>& mPeriodicIdVar;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /// Duplicate EquationIds to the second node on each periodic pair
    template< class TVariable >
    void CopyEquationId(const TVariable& rVariable,
                        ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
        rDest.pGetDof( rVariable )->SetEquationId( rOrigin.pGetDof(rVariable)->EquationId() );
    }

    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ResidualBasedBlockBuilderAndSolverPeriodic */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif // KRATOS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H
