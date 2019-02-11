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

#ifndef KRATOS_GENERAL_CONVERGENCE_CRITERIA_H
#define KRATOS_GENERAL_CONVERGENCE_CRITERIA_H

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
///@addtogroup IncompressibleFluidApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for fluid problems.
/**
 This class implements a convergence control based on nodal velocity and
 pressure values. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template <class TSparseSpace, class TDenseSpace>
class GeneralConvergenceCriteria : public ConvergenceCriteria<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GeneralConvergenceCriteria);

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef OpenMPUtils::PartitionVector PartitionVector;

    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param VelRatioTolerance Relative tolerance for velocity error
     * @param VelAbsTolerance Absolute tolerance for velocity error
     * @param PrsRatioTolerance Relative tolerance for presssure error
     * @param PrsAbsTolerance Absolute tolerance for presssure error
     */
    GeneralConvergenceCriteria(std::map<KeyType, TDataType>& rRatioTolerance,
                               std::map<KeyType, TDataType>& rAbsTolerance)
        : ConvergenceCriteria<TSparseSpace, TDenseSpace>(),
          mRatioTolerance(rRatioTolerance),
          mAbsTolerance(rAbsTolerance)
    {
    }

    /// Destructor.
    ~GeneralConvergenceCriteria() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Compute relative and absoute error.
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(ModelPart& rModelPart,
                      DofsArrayType& rDofSet,
                      const TSystemMatrixType& A,
                      const TSystemVectorType& Dx,
                      const TSystemVectorType& b) override
    {
        if (SparseSpaceType::Size(Dx) != 0) // if we are solving for something
        {
            int NumDofs = rDofSet.size();

            // Initialize
            std::map<KeyType, TDataType> solution_norm;
            std::map<KeyType, TDataType> increase_norm;
            std::map<KeyType, unsigned int> dof_num;

            for (auto variable_data : mRatioTolerance)
            {
                solution_norm[variable_data.first] = 0.0;
                dof_num[variable_data.first] = 0;
            }

            for (auto variable_data : mRatioTolerance)
                increase_norm[variable_data.first] = 0.0;

            // Set a partition for OpenMP
            PartitionVector DofPartition;
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(NumDofs, NumThreads, DofPartition);

// Loop over Dofs
#pragma omp parallel reduction(+:VelSolutionNorm,PrSolutionNorm,VelIncreaseNorm,PrIncreaseNorm,VelDofNum,PrDofNum)
            {
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin =
                    rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd =
                    rDofSet.begin() + DofPartition[k + 1];

                std::size_t DofId;
                TDataType DofValue;
                TDataType DofIncr;

                for (typename DofsArrayType::iterator itDof = DofBegin;
                     itDof != DofEnd; ++itDof)
                {
                    if (itDof->IsFree())
                    {
                        DofId = itDof->EquationId();
                        DofValue = itDof->GetSolutionStepValue(0);
                        DofIncr = Dx[DofId];

                        KeyType CurrVar = itDof->GetVariable().Key();

                        solution_norm[CurrVar] += DofValue * DofValue;
                        increase_norm[CurrVar] += DofIncr * DofIncr;
                        dof_num[CurrVar] += 1;
                    }
                }
            }

            for (auto variable_data : mRatioTolerance)
                if (solution_norm[variable_data.first] == 0.0)
                    solution_norm[variable_data.first] = 1.0;

            std::map<KeyType, TDataType> ratio;
            std::map<KeyType, TDataType> ratio_abs;
            for (auto variable_data : mRatioTolerance)
            {
                KeyType curr_var = variable_data.first;
                ratio[curr_var] = increase_norm[curr_var] / solution_norm[curr_var];
                ratio_abs[curr_var] = sqrt(increase_norm[curr_var]) /
                                      static_cast<int>(dof_num[curr_var]);
            }

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "CONVERGENCE CHECK:" << std::endl;
                for (auto variable_data : mRatioTolerance)
                {
                    KeyType cur_var = variable_data.first;
                    std::cout << "   " << cur_var;
                    std::cout << ": ratio = " << ratio[cur_var]
                              << "; exp.ratio = " << mRatioTolerance[cur_var]
                              << std::endl;
                    std::cout << "   " << cur_var;
                    std::cout << ": abs   = " << ratio_abs[cur_var]
                              << ";   exp.abs = " << mAbsTolerance[cur_var]
                              << std::endl;
                }
            }

            for (auto variable_data : mAbsTolerance)
            {
                KeyType cur_var = variable_data.first;
                if ((ratio[cur_var] > mRatioTolerance[cur_var]) &&
                    (ratio_abs[cur_var] > mAbsTolerance[cur_var]))
                    return false;
            }

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
            }

            return true;
        }
        else // in this case all the displacements are imposed!
        {
            return true;
        }
    }

    /// Initialize this class before using it
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem. (unused)
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                DofsArrayType& rDofSet,
                                const TSystemMatrixType& A,
                                const TSystemVectorType& Dx,
                                const TSystemVectorType& b) override
    {
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              DofsArrayType& rDofSet,
                              const TSystemMatrixType& A,
                              const TSystemVectorType& Dx,
                              const TSystemVectorType& b) override
    {
    }

    ///@} // Operations

private:
    std::map<KeyType, TDataType>& mRatioTolerance;
    std::map<KeyType, TDataType>& mAbsTolerance;
};

///@} // Kratos classes

///@} // Application group
} // namespace Kratos

#endif /* _VEL_PR_CRITERIA_H */
