//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Natalia Saiapova
//

#if !defined( KRATOS_GENERAL_CONVERGENCE_CRITERIA_H_INCLUDED )
#define KRATOS_GENERAL_CONVERGENCE_CRITERIA_H_INCLUDED

// System includes
#include <unordered_map>
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/color_utilities.h"
#include "input_output/logger_table_output.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
template<class TSparseSpace,
         class TDenseSpace>
class GeneralConvergenceCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeneralConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION(GeneralConvergenceCriteria);

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef VariableData::KeyType KeyType;

    typedef Variable<double> DoubleVariableType;

    typedef Variable< array_1d< double, 3> > Array3VariableType;

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;

    typedef std::size_t SizeType;

    typedef std::size_t IndexType;

    ///@}
    ///@name  Enum's
    ///@{

    enum class BasisVectorType
    {
        RESIDUAL,
        SOLUTION_UPDATE
    };

    ///@}
    ///@name Life Cycle
    ///@{

    ///
    GeneralConvergenceCriteria(
        TDataType RelativeTolerance,
        TDataType AbsoluteTolerance,
        const std::string& rBasisVectorType,
        const bool PrintColors = false)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mPrintColors(PrintColors)
    {
        SelectBasisVectorType(rBasisVectorType);

        mRelativeTolerances.resize(1);
        mAbsoluteTolerances.resize(1);

        mRelativeTolerances[0] = RelativeTolerance;
        mAbsoluteTolerances[0] = AbsoluteTolerance;
    }

    ///
    GeneralConvergenceCriteria(Parameters ThisParameters)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        const Parameters default_params( R"({
            "basis_vector_type"     : "residual",
            "variables_to_separate" : [],
            "relative_tolerances"   : [],
            "absolute_tolerances"   : [],
            "other_dofs_name"       : "Other-Dofs",
            "print_colors"          : false,
            "echo_level"            : 0
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_params);

        mPrintColors = ThisParameters["print_colors"].GetBool();
        SelectBasisVectorType(ThisParameters["basis_vector_type"].GetString());
        this->SetEchoLevel(ThisParameters["echo_level"].GetInt());

        const SizeType num_vars_to_separate = ThisParameters["variables_to_separate"].size() + 1; // +1 bcs of the other dofs

        // Size check // TODO add "+1" in description!
        KRATOS_ERROR_IF(num_vars_to_separate != ThisParameters["relative_tolerances"].size())
            << "Your list of variables is not the same size as the list of "
            << "relative_tolerances" << std::endl;
        KRATOS_ERROR_IF(num_vars_to_separate != ThisParameters["absolute_tolerances"].size())
            << "Your list of variables is not the same size as the list of "
            << "absolute_tolerances" << std::endl;

        mRelativeTolerances.resize(num_vars_to_separate);
        mAbsoluteTolerances.resize(num_vars_to_separate);

        mVariableNames.resize(num_vars_to_separate);
        mVariableNames[num_vars_to_separate-1] = ThisParameters["other_dofs_name"].GetString();

        for (IndexType i_var=0; i_var<num_vars_to_separate; ++i_var) {
            mRelativeTolerances[i_var] = ThisParameters["relative_tolerances"].GetArrayItem(i_var).GetDouble();
            mAbsoluteTolerances[i_var] = ThisParameters["absolute_tolerances"].GetArrayItem(i_var).GetDouble();
        }

        for (IndexType i_var=0; i_var<num_vars_to_separate-1; ++i_var) {
            const std::string& r_variable_name = ThisParameters["variables_to_separate"].GetArrayItem(i_var).GetString();
            KeyType variable_key;

            mVariableNames[i_var] = r_variable_name;

            if (KratosComponents<DoubleVariableType>::Has(r_variable_name)) {
                variable_key = KratosComponents< DoubleVariableType >::Get(r_variable_name).Key();
                mKeyToIndexMap[variable_key] = i_var;
            } else if (KratosComponents< Array3VariableType >::Has(r_variable_name)) {
                // In this case all the variables point to the same index
                variable_key = KratosComponents< ComponentVariableType >::Get(r_variable_name+std::string("_X")).Key();
                mKeyToIndexMap[variable_key] = i_var;
                variable_key = KratosComponents< ComponentVariableType >::Get(r_variable_name+std::string("_Y")).Key();
                mKeyToIndexMap[variable_key] = i_var;
                variable_key = KratosComponents< ComponentVariableType >::Get(r_variable_name+std::string("_Z")).Key();
                mKeyToIndexMap[variable_key] = i_var;
            }
            else if (KratosComponents< ComponentVariableType >::Has(r_variable_name)) {
                variable_key = KratosComponents< ComponentVariableType >::Get(r_variable_name).Key();
                mKeyToIndexMap[variable_key] = i_var;
            }
            else {
                KRATOS_ERROR << "Only Double (e.g. PRESSURE), Array3D (e.g. VELOCITY) or Component "
                    << "(e.g. DISPLACEMENT_X) variables are allowed in the variables list" << std::endl;
            }
        }

        // TODO test if a component and its Array are specified!
    }

    /// Destructor.
    ~GeneralConvergenceCriteria() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
    }

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
        const TSystemVectorType& r_vec = (mBasisVectorType == BasisVectorType::RESIDUAL) ? b : Dx;

        // first check if we are solving for something
        if (SparseSpaceType::Size(r_vec) == 0) return true;

        const SizeType num_vars_to_separate = mVariableNames.size();

        // Initialize the vectors holding the residuals
        std::vector<TDataType> relative_residuals(num_vars_to_separate, TDataType());
        std::vector<TDataType> absolute_residuals(num_vars_to_separate, TDataType());
        std::vector<int> num_dofs(num_vars_to_separate, int());

        // Calculate the residuals
        const bool is_mpi_execution = rModelPart.GetCommunicator().TotalProcesses() > 1;
        const bool separate_variables = num_vars_to_separate > 1;

        if (is_mpi_execution) {
            const double rank = rModelPart.GetCommunicator().MyPID(); // To compare with PARTITION_INDEX, which is a double variable
            if (separate_variables) {
                CalculateSeparatedResidualsMPI(rDofSet, r_vec, relative_residuals, absolute_residuals, num_dofs, rank);
            } else {
                CalculateResidualsMPI(rDofSet, r_vec, relative_residuals, absolute_residuals, num_dofs, rank);
            }
        } else {
            if (separate_variables) {
                CalculateSeparatedResiduals(rDofSet, r_vec, relative_residuals, absolute_residuals, num_dofs);
            } else {
                CalculateResiduals(rDofSet, r_vec, relative_residuals, absolute_residuals, num_dofs);
            }
        }

        FinalizeResidualsCalculation(relative_residuals, absolute_residuals, num_dofs);

        // This is done on each rank, since the residuals were synchronized before
        std::vector<bool> conv_vec(num_vars_to_separate); // save the convegence info for plotting
        bool is_converged = true;

        for (SizeType i=0; i<num_vars_to_separate; ++i) {
            // Check each variable for convergence
            conv_vec[i] = relative_residuals[i] < mRelativeTolerances[i] ||
                          absolute_residuals[i] < mAbsoluteTolerances[i];

            is_converged = is_converged && conv_vec[i]; // Check overall convergence
        }

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            const int nonlin_iteration_number = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            PrintConvergenceInfo(is_converged, conv_vec, relative_residuals, absolute_residuals, nonlin_iteration_number);
        }

        return is_converged;
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart
     * @return 0 all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Check(rModelPart);

        if (rModelPart.GetCommunicator().TotalProcesses() > 1) { // mpi-execution
            KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX))
                << "PARTITION_INDEX is not a solutionstep-variable!" << std::endl;
        }

        // TODO check if the variables that are being checked here are in the ModelPart!
        // TODO check if the variables that are being checked here Dofs!

        return 0;
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "GeneralConvergenceCriteria";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GeneralConvergenceCriteria";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<TDataType> mRelativeTolerances;
    std::vector<TDataType> mAbsoluteTolerances;

    std::unordered_map<KeyType, IndexType> mKeyToIndexMap;
    std::vector<std::string> mVariableNames;

    BasisVectorType mBasisVectorType;

    bool mPrintColors = false;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateResiduals(
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rVector,
        std::vector<TDataType>& rRelativeResiduals,
        std::vector<TDataType>& rAbsoluteResiduals,
        std::vector<int>& rNumDofs)
    {
        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(rDofSet.size()); ++i) {
            TDataType dof_value;
            TDataType dof_incr;

            const auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree()) {
                dof_value = it_dof->GetSolutionStepValue(0);
                dof_incr = TSparseSpace::GetValue(rVector, it_dof->EquationId());

                #pragma omp critical
                {
                    rRelativeResiduals[0] += dof_value * dof_value;
                    rAbsoluteResiduals[0] += dof_incr * dof_incr;
                    ++rNumDofs[0];
                }
            }
        }
    }

    void CalculateSeparatedResiduals(
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rVector,
        std::vector<TDataType>& rRelativeResiduals,
        std::vector<TDataType>& rAbsoluteResiduals,
        std::vector<int>& rNumDofs)
    {
        const IndexType other_dofs_vec_index = mVariableNames.size()-1;

        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(rDofSet.size()); ++i) {
            IndexType vec_index;
            TDataType dof_value;
            TDataType dof_incr;

            const auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree()) {
                dof_value = it_dof->GetSolutionStepValue(0);
                dof_incr = TSparseSpace::GetValue(rVector, it_dof->EquationId());

                KeyType dof_var_key = it_dof->GetVariable().Key();

                // Here we are getting the index that belongs to the corresponding variable key
                // If the key for which we want to get the index does not exist, we get 0,
                // which corresponds to the "remaining" dofs
                // at and count have constant (worst case linear) complexity, so this should be fine since the map is small
                vec_index = (mKeyToIndexMap.count(dof_var_key)) ? mKeyToIndexMap.at(dof_var_key) : other_dofs_vec_index;

                #pragma omp critical
                {
                    rRelativeResiduals[vec_index] += dof_value * dof_value;
                    rAbsoluteResiduals[vec_index] += dof_incr * dof_incr;
                    ++rNumDofs[vec_index];
                }
            }
        }
    }

    void CalculateResidualsMPI(
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rVector,
        std::vector<TDataType>& rRelativeResiduals,
        std::vector<TDataType>& rAbsoluteResiduals,
        std::vector<int>& rNumDofs,
        const double Rank)
    {
        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(rDofSet.size()); ++i) {
            TDataType dof_value;
            TDataType dof_incr;

            const auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX)==Rank) {
                dof_value = it_dof->GetSolutionStepValue(0);
                dof_incr = TSparseSpace::GetValue(rVector, it_dof->EquationId());

                #pragma omp critical
                {
                    rRelativeResiduals[0] += dof_value * dof_value;
                    rAbsoluteResiduals[0] += dof_incr * dof_incr;
                    ++rNumDofs[0];
                }
            }
        }
        SynchronizeResiduals(rRelativeResiduals, rAbsoluteResiduals, rNumDofs);
    }

    void CalculateSeparatedResidualsMPI(
        const DofsArrayType& rDofSet,
        const TSystemVectorType& rVector,
        std::vector<TDataType>& rRelativeResiduals,
        std::vector<TDataType>& rAbsoluteResiduals,
        std::vector<int>& rNumDofs,
        const double Rank)
    {
        const IndexType other_dofs_vec_index = mVariableNames.size()-1;

        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(rDofSet.size()); ++i) {
            IndexType vec_index;
            TDataType dof_value;
            TDataType dof_incr;

            const auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX)==Rank) {
                dof_value = it_dof->GetSolutionStepValue(0);
                dof_incr = TSparseSpace::GetValue(rVector, it_dof->EquationId());

                KeyType dof_var_key = it_dof->GetVariable().Key();

                // Here we are getting the index that belongs to the corresponding variable key
                // If the key for which we want to get the index does not exist, we get 0,
                // which corresponds to the "remaining" dofs
                // at and count have constant (worst case linear) complexity, so this should be fine since the map is small
                vec_index = (mKeyToIndexMap.count(dof_var_key)) ? mKeyToIndexMap.at(dof_var_key) : other_dofs_vec_index;

                #pragma omp critical
                {
                    rRelativeResiduals[vec_index] += dof_value * dof_value;
                    rAbsoluteResiduals[vec_index] += dof_incr * dof_incr;
                    ++rNumDofs[vec_index];
                }
            }
        }
        SynchronizeResiduals(rRelativeResiduals, rAbsoluteResiduals, rNumDofs);
    }

    void SynchronizeResiduals(
        std::vector<TDataType>& rRelativeResiduals,
        std::vector<TDataType>& rAbsoluteResiduals,
        std::vector<int>& rNumDofs)
    {
        const SizeType num_vars_to_separate = rNumDofs.size();

        // concatenate the vectors to have only one call to MPI
        std::vector<TDataType> residuals = rRelativeResiduals;
        residuals.reserve(2 * num_vars_to_separate);
        std::copy(rAbsoluteResiduals.begin(), rAbsoluteResiduals.end(), residuals.begin() + num_vars_to_separate);

        // Synchroizing them across ranks // TODO use data-communicator
        // rModelPart.GetCommunicator().SumAll(residuals);
        // rModelPart.GetCommunicator().SumAll(rNumDofs);

        // Then afterwards split them again
        std::copy(residuals.begin(), residuals.begin() + num_vars_to_separate, rRelativeResiduals.begin());
        std::copy(residuals.begin() + num_vars_to_separate, residuals.end(), rAbsoluteResiduals.begin());
    }

    void FinalizeResidualsCalculation(
        std::vector<TDataType>& rRelativeResiduals,
        std::vector<TDataType>& rAbsoluteResiduals,
        std::vector<int>& rNumDofs)
    {
        auto sqrt_fct = [](const TDataType & el) -> TDataType { return std::sqrt(el); };
        // Finish applying the L2-Norm
        std::transform(rRelativeResiduals.begin(), rRelativeResiduals.end(), rRelativeResiduals.begin(), sqrt_fct);
        std::transform(rAbsoluteResiduals.begin(), rAbsoluteResiduals.end(), rAbsoluteResiduals.begin(), sqrt_fct);
        // Take into account the size of the domain
        std::transform(rNumDofs.begin(), rNumDofs.end(), rNumDofs.begin(), sqrt_fct);

        const SizeType num_vars_to_separate = rNumDofs.size();

        for (SizeType i=0; i<num_vars_to_separate; ++i) {
            if (IsAlmostZero(rRelativeResiduals[i])) {
                rRelativeResiduals[i] = 1.0;
            }

            if (rNumDofs[i] == 0) {
                rNumDofs[i] = 1; // this might be the case if no "other" dofs are present
            }

            // Compute the final residuals
            rRelativeResiduals[i]  = rAbsoluteResiduals[i] / rRelativeResiduals[i];
            rAbsoluteResiduals[i] /= rNumDofs[i];
        }
    }

    bool IsAlmostZero(
        const double Value,
        const double Eps = std::numeric_limits<double>::epsilon()) const
    {
        return std::abs(Value) <= Eps * std::abs(Value);
    }

    // Function to select the vector to be checked for convergence
    void SelectBasisVectorType(const std::string& rBasisVectorType)
    {
        if (rBasisVectorType == "residual") {
            mBasisVectorType = BasisVectorType::RESIDUAL;
        } else if (rBasisVectorType == "solution_update") {
            mBasisVectorType = BasisVectorType::SOLUTION_UPDATE;
        }
        else {
            KRATOS_ERROR << "Wrong \"BasisVectorType\", use \"residual\" "
                << "or \"solution_update\"" << std::endl;
        }
    }

    /*
    Function to print info abt the convergence
    Note that we are flushing the buffer at the end only!
    */
    void PrintConvergenceInfo(const bool IsConverged,
                              const std::vector<bool>& rConvergenceInfoVector,
                              const std::vector<TDataType>& rRelativeResiduals,
                              const std::vector<TDataType>& rAbsoluteResiduals,
                              const int NonlinIterationNumber) const
    {
        KRATOS_INFO("ConvergenceCriteria") << "Convergence Check; Iteration "
            << NonlinIterationNumber << "\n";

        for (IndexType i=0; i<mVariableNames.size(); ++i) {
            std::stringstream conv_info;
            if (rConvergenceInfoVector[i]) {
                if (mPrintColors) conv_info << BOLDFONT(FGRN("converged"));
                else conv_info << "converged";
            } else {
                if (mPrintColors) conv_info << BOLDFONT(FRED("not converged"));
                else conv_info << "not converged";
            }

            KRATOS_INFO("ConvergenceCriteria") << "\t" << mVariableNames[i] << ": "<< conv_info.str()
                << " | ratio = " << rRelativeResiduals[i] << "; exp.ratio = " << mRelativeTolerances[i] << " | "
                << "abs = "      << rAbsoluteResiduals[i] << "; exp.abs = "   << mAbsoluteTolerances[i] << "\n";
        }

        if (IsConverged) {
            if (mPrintColors) {
                KRATOS_INFO("ConvergenceCriteria") << BOLDFONT(FGRN("Convergence is achieved in Iteration "))
                    << NonlinIterationNumber << std::endl;
            } else {
                KRATOS_INFO("ConvergenceCriteria") << "Convergence is achieved in Iteration "
                    << NonlinIterationNumber << std::endl;
            }
        } else {
            if (mPrintColors) {
            KRATOS_INFO("ConvergenceCriteria") << BOLDFONT(FRED("Convergence is not achieved in Iteration "))
                << NonlinIterationNumber << std::endl;
            } else {
            KRATOS_INFO("ConvergenceCriteria") << "Convergence is not achieved in Iteration "
                << NonlinIterationNumber << std::endl;
            }
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GeneralConvergenceCriteria& operator=(GeneralConvergenceCriteria const& rOther){}

    /// Copy constructor.
    GeneralConvergenceCriteria(GeneralConvergenceCriteria const& rOther){}

    ///@}

}; // Class GeneralConvergenceCriteria

///@}
///@name Input and output
///@{

// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 GeneralConvergenceCriteria& rThis){}

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const GeneralConvergenceCriteria& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GENERAL_CONVERGENCE_CRITERIA_H_INCLUDED  defined


