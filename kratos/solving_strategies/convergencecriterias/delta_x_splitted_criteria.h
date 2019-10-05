//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on the work of Jordi Cotela and Vicente Mataix Ferrandiz)
//

#if !defined(KRATOS_DELTA_X_SPLITTED_CRITERIA_H_INCLUDED )
#define  KRATOS_DELTA_X_SPLITTED_CRITERIA_H_INCLUDED


// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class DeltaXSplittedCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DeltaXSplittedCriteria
    KRATOS_CLASS_POINTER_DEFINITION(DeltaXSplittedCriteria);

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename DofsArrayType::TDataType DofType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DeltaXSplittedCriteria(
        const std::string& rMainDofsName,
        const std::string& rSplittedDofsName,
        const std::unordered_set<std::string>& rMainDofNames, // using set to avoid duplications
        const TDataType MainDofsAbsTolerance,
        const TDataType MainDofsRatioTolerance,
        const TDataType SplittedDofsAbsTolerance,
        const TDataType SplittedDofsRatioTolerance) :
            mMainDofsName(rMainDofsName),
            mSplittedDofsName(rSplittedDofsName),
            mMainDofsAbsTolerance(MainDofsAbsTolerance),
            mMainDofsRatioTolerance(MainDofsRatioTolerance),
            mSplittedDofsAbsTolerance(SplittedDofsAbsTolerance),
            mSplittedDofsRatioTolerance(SplittedDofsRatioTolerance)
    {
        mMainDofVariables.clear();
        for (const auto& rDofVarName : rMainDofNames) {
            mMainDofVariables.push_back(&KratosComponents<VariableData>::Get(rDofVarName));
        }
    }

    /// Destructor.
    virtual ~DeltaXSplittedCriteria(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
        mRank = rModelPart.GetCommunicator().MyPID();
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @return 0 all OK, 1 otherwise
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(TSparseSpace::IsDistributed() == rModelPart.IsDistributed()) << "Mismatch between ModelPart and ConvergenceCriteria IsDistributed!" << std::endl;
        if (TSparseSpace::IsDistributed()) {
            KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX)) << "PARTITION_INDEX nneds to be a solution-step variable on Modelpart: " << rModelPart.Name();
        }

        return BaseType::Check(rModelPart);
        KRATOS_CATCH("");
    }

    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        if (rDx.size() == 0) { // check if we are solving for something
            return;
        }

        std::vector<TDataType> norms(0.0, 4);
        std::vector<SizeType> num_dofs(0, 2);

        ComputeNorms(rDofSet, rDx, norms, num_dofs);

        std::vector<TDataType> global_norms(0.0, 4);
        std::vector<SizeType> global_num_dofs(0, 2);

        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();

        r_data_comm.Sum(norms, global_norms, 0);
        r_data_comm.Sum(num_dofs, global_num_dofs, 0);

        int is_converged(0);

        if (mRank == 0) {
            // if(VelSolutionNorm == 0.0)
            //     VelSolutionNorm = 1.0;
            // if(PrSolutionNorm == 0.0)
            //     PrSolutionNorm = 1.0;

            // ratio_main     = std::sqrt(VelIncreaseNorm/VelSolutionNorm);
            // ratio_splitted = std::sqrt(PrIncreaseNorm/PrSolutionNorm);

            // abs_main     = std::sqrt(VelIncreaseNorm)/ static_cast<TDataType>(global_num_dofs[0]);
            // abs_splitted = std::sqrt(PrIncreaseNorm)/ static_cast<TDataType>(global_num_dofs[1]);

            // is_converged = (ratio_main   <= mMainDofsRatioTolerance ||
            //                 abs_main     <= mMainDofsAbsTolerance) &&
            //                (ratio_splitted <= mSplittedDofsRatioTolerance ||
            //                 abs_splitted   <= mSplittedDofsAbsTolerance)
            KRATOS_INFO_IF("DeltaXSplittedCriteria", is_converged && this->GetEchoLevel() > 0) << "Convergence is achieved" << std::endl;
        }

        // scatter is_converged

        r_data_comm.Broadcast(is_converged, 0);

        return static_cast<bool>(is_converged);
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DeltaXSplittedCriteria" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "DeltaXSplittedCriteria";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

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


    ///@}
    ///@name Member Variables
    ///@{

    std::string mMainDofsName;
    std::string mSplittedDofsName;
    TDataType mMainDofsAbsTolerance;
    TDataType mMainDofsRatioTolerance;
    TDataType mSplittedDofsAbsTolerance;
    TDataType mSplittedDofsRatioTolerance;
    std::vector<VariableData*> mMainDofVariables;
    int mRank = 0;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    inline bool ConsiderDof(const DofType& rDof) const
    {
        if (TSparseSpace::IsDistributed()) {
            return rDof.IsFree() && rDof.GetSolutionStepValue(PARTITION_INDEX) == mRank;
        } else {
            return rDof.IsFree();
        }
    }

    void ComputeNorms(
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rDx,
        std::vector<TDataType>& rNorms,
        std::vector<TDataType>& rNumDofs) const
    {
        TDataType increase_norm_main_dofs(0.0);
        TDataType solution_norm_main_dofs(0.0);
        TDataType increase_norm_splitted_dofs(0.0);
        TDataType solution_norm_splitted_dofs(0.0);
        SizeType num_main_dofs(0);
        SizeType num_splitted_dofs(0);

        // Loop over Dofs
        #pragma omp parallel for reduction(+:increase_norm_main_dofs, solution_norm_main_dofs, increase_norm_splitted_dofs, solution_norm_splitted_dofs, num_main_dofs, num_splitted_dofs)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
            auto it_dof = rDofSet.begin() + i;

            if (ConsiderDof(it_dof)) {
                const TDataType dof_incr = TSparseSpace::GetValue(rDx, it_dof->EquationId());
                const TDataType dof_val  = it_dof->GetSolutionStepValue(0);

                if (std::find(mMainDofVariables.begin(), mMainDofVariables.end(), it_dof->GetVariable())) {
                    num_main_dofs++;
                    increase_norm_main_dofs += std::pow(dof_incr, 2);
                    solution_norm_main_dofs += std::pow(dof_val, 2);
                } else {
                    num_splitted_dofs++;
                    increase_norm_splitted_dofs += std::pow(dof_incr, 2);
                    solution_norm_splitted_dofs += std::pow(dof_val, 2);
                }
            }
        }
        rNorms[0] = increase_norm_main_dofs;
        rNorms[1] = solution_norm_main_dofs;
        rNorms[2] = increase_norm_splitted_dofs;
        rNorms[3] = solution_norm_splitted_dofs;

        rNumDofs[0] = num_main_dofs;
        rNumDofs[1] = num_splitted_dofs;
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DeltaXSplittedCriteria& operator=(DeltaXSplittedCriteria const& rOther){}

    /// Copy constructor.
    DeltaXSplittedCriteria(DeltaXSplittedCriteria const& rOther){}

    ///@}

}; // Class DeltaXSplittedCriteria

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_DELTA_X_SPLITTED_CRITERIA_H_INCLUDED defined


