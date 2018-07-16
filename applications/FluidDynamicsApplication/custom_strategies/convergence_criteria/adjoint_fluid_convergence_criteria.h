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


#ifndef KRATOS_ADJOINT_FLUID_CONV_CRITERIA_H
#define	KRATOS_ADJOINT_FLUID_CONV_CRITERIA_H

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@addtogroup IncompressibleFluidApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for fluid problems.
/**
 This class implements a convergence control based on nodal adjoint velocity and
 adjoint pressure values. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template<   class TSparseSpace,
            class TDenseSpace >
class AdjointFluidConvergenceCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( AdjointFluidConvergenceCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

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
    AdjointFluidConvergenceCriteria(  TDataType AdjointVariableRatioTolerance,
                                      TDataType AdjointVariableAbsTolerance)
            : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mAdjointVariableRatioTolerance = AdjointVariableRatioTolerance;
        mAdjointVariableAbsTolerance = AdjointVariableAbsTolerance;

    }

    /// Destructor.
    ~AdjointFluidConvergenceCriteria() override {}

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
    bool PostCriteria(  ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        const TSystemMatrixType& A,
                        const TSystemVectorType& Dx,
                        const TSystemVectorType& b ) override
    {
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            // Initialize
            TDataType AdjointVariableSolutionNorm = 0.0;
            TDataType AdjointVariableIncreaseNorm = 0.0;
            unsigned int AdjointVariableDofNum(0);

            // Set a partition for OpenMP
            int NumDofs = rDofSet.size();
            PartitionVector DofPartition;
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(NumDofs,NumThreads,DofPartition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:AdjointVariableSolutionNorm,AdjointVariableIncreaseNorm,AdjointVariableDofNum)
            {
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k+1];

                std::size_t DofId;
                TDataType DofValue;
                TDataType DofIncr;

                for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
                {
                    if (itDof->IsFree())
                    {
                        DofId = itDof->EquationId();
                        DofValue = itDof->GetSolutionStepValue(0);
                        DofIncr = Dx[DofId];

                        AdjointVariableSolutionNorm += DofValue * DofValue;
                        AdjointVariableIncreaseNorm += DofIncr * DofIncr;
                        ++AdjointVariableDofNum;
                    }
                }
            }

            if(AdjointVariableDofNum == 0.0)
                AdjointVariableDofNum = 1.0;

            TDataType AdjointVariableRatio = sqrt(AdjointVariableIncreaseNorm/AdjointVariableSolutionNorm);
            TDataType AdjointVariableAbs = sqrt(AdjointVariableIncreaseNorm)/ static_cast<TDataType>(AdjointVariableDofNum);

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                int iteration = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                std::cout << "Non linear Iteration: "<<iteration<<"; Adjoint variable convergence: Ratio = " << AdjointVariableRatio <<"; exp.ratio = " << mAdjointVariableRatioTolerance << " abs = " << AdjointVariableAbs << " exp.abs = " << mAdjointVariableAbsTolerance << std::endl;
            }

            if (    (AdjointVariableRatio <= mAdjointVariableRatioTolerance || AdjointVariableAbs <= mAdjointVariableAbsTolerance) )            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    std::cout << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        else //in this case all the displacements are imposed!
        {
            return true;
        }
    }

    /// Initialize this class before using it
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem. (unused)
     */
    void Initialize( ModelPart& rModelPart	) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(    ModelPart& rModelPart,
                                    DofsArrayType& rDofSet,
                                    const TSystemMatrixType& A,
                                    const TSystemVectorType& Dx,
                                    const TSystemVectorType& b ) override
    {}

    void FinalizeSolutionStep(  ModelPart& rModelPart,
                                DofsArrayType& rDofSet,
                                const TSystemMatrixType& A,
                                const TSystemVectorType& Dx,
                                const TSystemVectorType& b ) override
    {}

    ///@} // Operations

private:

    TDataType mAdjointVariableRatioTolerance;
    TDataType mAdjointVariableAbsTolerance;
};

///@} // Kratos classes

///@} // Application group
}

#endif	/* _VEL_PR_CRITERIA_H */
