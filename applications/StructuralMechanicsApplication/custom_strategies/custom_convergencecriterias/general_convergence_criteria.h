//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Natalia Saiapova
//                   Philipp Bucher
//


#if !defined(KRATOS_GENERAL_CONVERGENCE_CRITERIA_H_INCLUDED )
#define  KRATOS_GENERAL_CONVERGENCE_CRITERIA_H_INCLUDED

/*
Here is a list of files that you can look at for reference:
- applications/StructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/displacement_and_other_dof_criteria.h
- applications/StructuralMechanicsApplication/custom_strategies/custom_convergencecriterias/residual_displacement_and_other_dof_criteria.h
- applications/FluidDynamicsApplication/custom_strategies/convergence_criteria/vel_pr_criteria.h
- applications/trilinos_application/custom_strategies/convergencecriterias/trilinos_displacement_criteria.h
- applications/trilinos_application/custom_strategies/convergencecriterias/trilinos_up_criteria.h
- kratos/solving_strategies/convergencecriterias/displacement_criteria.h
- kratos/solving_strategies/convergencecriterias/residual_criteria.h


*/

// System includes
#include <unordered_map>


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

    enum BasisType
    {
        RESIDUAL,
        SOLUTION_UPDATE
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeneralConvergenceCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm,
        Parameters rParameters = Parameters(R"({})"),
        const std::string OtherDofsName="Other Dofs")
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        Parameters default_params( R"({
            "basis_vector_type" : "residual",
            "variables_to_separate" : [],
            "relative_convergence_tolerances" : [],
            "absolut_convergence_tolerances" : []
        })" );

        rParameters.ValidateAndAssignDefaults(default_params);

        const std::string& r_basis_type = rParameters["basis_vector_type"].GetString();
        if (r_basis_type == "residual")             mBasisType = BasisType::RESIDUAL;
        else if (r_basis_type == "solution_update") mBasisType = BasisType::SOLUTION_UPDATE;
        else KRATOS_ERROR << "Wrong BasisType, use \"residual\" or \"solution_update\"" << std::endl;

        const SizeType n_variables = rParameters["variables_to_separate"].size();
        const SizeType num_vars_to_separate = n_variables + 1; // +1 bcs the "remaining" dofs are at pos 0
        const SizeType n_rel_conv_tol = rParameters["relative_convergence_tolerances"].size();
        const SizeType n_abs_conv_tol = rParameters["absolut_convergence_tolerances"].size();

        // Size check
        KRATOS_ERROR_IF(n_variables != n_rel_conv_tol) << "Your list of variables is not the same size as the list of relative_convergence_tolerances" << std::endl;
        KRATOS_ERROR_IF(n_variables != n_abs_conv_tol) << "Your list of variables is not the same size as the list of absolut_convergence_tolerances" << std::endl;

        mRatioTolerances.resize(num_vars_to_separate);
        mAbsTolerances.resize(num_vars_to_separate);

        mRatioResiduals.resize(num_vars_to_separate);
        mAbsResiduals.resize(num_vars_to_separate);
        mNumDofs.resize(num_vars_to_separate);

        mIndexToNameMap[0] = (num_vars_to_separate > 1) ? OtherDofsName + std::string(": ") : ""; // Print nothing in case of no separation

        mRatioTolerances[0] = NewRatioTolerance;
        mAbsTolerances[0] = AlwaysConvergedNorm;

        for (IndexType i_var = 1; i_var < num_vars_to_separate; ++i_var)
        {
            const IndexType reverse_idx = num_vars_to_separate - i_var - 1; // For printing in proper order
            mRatioTolerances[i_var] = rParameters["relative_convergence_tolerances"].GetArrayItem(reverse_idx).GetDouble();
            mAbsTolerances[i_var] = rParameters["absolut_convergence_tolerances"].GetArrayItem(reverse_idx).GetDouble();

            const std::string& variable_name = rParameters["variables_to_separate"].GetArrayItem(reverse_idx).GetString();
            KeyType the_key;

            mIndexToNameMap[i_var] = variable_name + std::string(": ") ;

            if (KratosComponents<DoubleVariableType>::Has(variable_name))
            {
                the_key = KratosComponents< DoubleVariableType >::Get(variable_name).Key();
                mKeyToIndexMap[the_key] = i_var;
            }
            else if (KratosComponents< Array3VariableType >::Has(variable_name))
            {
                // In this case all the variables point to the same index
                the_key = KratosComponents< ComponentVariableType >::Get(variable_name+std::string("_X")).Key();
                mKeyToIndexMap[the_key] = i_var;
                the_key = KratosComponents< ComponentVariableType >::Get(variable_name+std::string("_Y")).Key();
                mKeyToIndexMap[the_key] = i_var;
                the_key = KratosComponents< ComponentVariableType >::Get(variable_name+std::string("_Z")).Key();
                mKeyToIndexMap[the_key] = i_var;
            }
            else if (KratosComponents< ComponentVariableType >::Has(variable_name)) //case of component variable)
            {
                the_key = KratosComponents< ComponentVariableType >::Get(variable_name).Key();
                mKeyToIndexMap[the_key] = i_var;
            }
            else
            {
                KRATOS_ERROR << "Only Double (e.g. PRESSURE), Array3D (e.g. VELOCITY) or Component "
                    << "(e.g. VELOCITY_X) variables are allowed in the variables list" << std::endl;
            }
        }

        // TODO test if a component and its Array are specified!


    }

    /// Destructor.
    virtual ~GeneralConvergenceCriteria(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //Criterias that need to be called after getting the solution
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
        // Natasha we will use this to print whether convergence has been achieved or not
        // I think it is quite cool :D
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << "Some regular print" << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << BOLDFONT("Some regular print") << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << FGRN("Some regular print") << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << BOLDFONT(FRED("   Not achieved")) << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << BOLDFONT(FBLU("   Not achieved")) << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << UNDL(FYEL("   Not achieved")) << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << UNDL(FMAG("   Not achieved")) << std::endl;
        // KRATOS_INFO("CONVERGENCE_CRITERIA") << UNDL(FCYN("   Not achieved")) << std::endl;

        const TSystemVectorType& r_vec = (mBasisType == BasisType::RESIDUAL) ? b : Dx;

        if (SparseSpaceType::Size(r_vec) != 0) // if we are solving for something
        {
            // Initialize the vectors
            std::fill(mRatioResiduals.begin(), mRatioResiduals.end(), TDataType());
            std::fill(mAbsResiduals.begin(), mAbsResiduals.end(), TDataType());
            std::fill(mNumDofs.begin(), mNumDofs.end(), SizeType());

            #pragma omp parallel for // comment bcs some vars have to be private to the thread
            for (int i = 0; i < static_cast<int>(rDofSet.size()); ++i)
            {
                IndexType dof_id;
                IndexType vec_index;
                TDataType dof_value;
                TDataType dof_incr;

                auto it_dof = rDofSet.begin() + i;

                if (it_dof->IsFree())
                {
                    dof_id = it_dof->EquationId();
                    dof_value = it_dof->GetSolutionStepValue(0);
                    dof_incr = r_vec[dof_id]; // TODO check if this is local => most probably implement a fct in the Space

                    KeyType dof_key = it_dof->GetVariable().Key();

                    // Here we are getting the index that belongs to the corresponding variable key
                    // If the key for which we want to get the index does not exist, we get 0,
                    // which corresponds to the "remaining" dofs
                    // at and count have constant (worst case linear) complexity, so this should be fine since the map is small
                    vec_index = (mKeyToIndexMap.count(dof_key)) ? mKeyToIndexMap.at(dof_key) : 0;

                    // These have to be atomic for omp if we don't do anything else
                    // Theoretically other operations are also possible
                    #pragma omp critical
                    {
                        mRatioResiduals[vec_index] += dof_value * dof_value;
                        mAbsResiduals[vec_index] += dof_incr * dof_incr;
                        ++mNumDofs[vec_index];
                    }
                }
            }

            // Computing the residuals

            const SizeType num_vars_to_separate = mRatioResiduals.size();

            // concatenate the vectors to have only one call to MPI
            std::vector<TDataType> residuals = mRatioResiduals;
            residuals.reserve(2 * num_vars_to_separate);
            std::copy(mAbsResiduals.begin(), mAbsResiduals.end(), residuals.begin() + num_vars_to_separate);

            // Synchroizing them across ranks
            rModelPart.GetCommunicator().SumAll(residuals);

            // Then afterwards split them again
            std::copy(residuals.begin(), residuals.begin() + num_vars_to_separate, mRatioResiduals.begin());
            std::copy(residuals.begin() + num_vars_to_separate, residuals.end(), mAbsResiduals.begin());

            auto sqrtFunction = [](const TDataType & el) -> TDataType { return std::sqrt(el); };
            // Finish applying the L2-Norm
            std::transform(mRatioResiduals.begin(), mRatioResiduals.end(), mRatioResiduals.begin(), sqrtFunction);
            std::transform(mAbsResiduals.begin(), mAbsResiduals.end(), mAbsResiduals.begin(), sqrtFunction);
            // Take into account the size of the domain
            std::transform(mNumDofs.begin(), mNumDofs.end(), mNumDofs.begin(), sqrtFunction);


            // Computing the final residuals and checking for convergence
            // This is done on each rank, since the residuals were synchronized before
            std::vector<bool> conv_vec(num_vars_to_separate); // save the convegence info for plotting
            bool is_converged = true;
            // @Natasha I changed this to a simple size-based loop, I think this is ok
            for (SizeType i=0; i<num_vars_to_separate; ++i)
            {
                // TODO implement isZero()
                if (mRatioResiduals[i] == 0.0) mRatioResiduals[i] = 1.0;
                if (mNumDofs[i] == 0) mNumDofs[i] = 0; // this might be the case if no "other" dofs are present

                // Compute the final residuals
                mRatioResiduals[i] = mAbsResiduals[i] / mRatioResiduals[i];
                mAbsResiduals[i] /= mNumDofs[i];

                // Check each variable for convergence
                conv_vec[i] = mRatioResiduals[i] < mRatioTolerances[i] ||
                              mAbsResiduals[i] < mAbsTolerances[i];

                is_converged = is_converged && conv_vec[i]; // Check overall convergence
            }

            // Printing information abt current residuals
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                const int nonlin_iteration_number = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                PrintConvergenceInfo(is_converged, conv_vec, nonlin_iteration_number);
            }

            return is_converged;
        }
        else
        {
            return true;
        }
    }

    void Initialize(
        ModelPart& rModelPart
    ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;

        // static std::stringstream buffer;
        // LoggerOutput::Pointer p_output(new LoggerTableOutput(buffer, {"Time Step", "Iteration Number        ", "Convergence        ", "Is converged"}));
        // Logger::AddOutput(p_output);
        // p_output->WriteHeader();

        // for (int i=0; i<5; ++i)
        // {

        //     LoggerMessage message1("Time Step");
        //     message1 << "445";
        //     LoggerMessage message2("Iteration Number");
        //     message2 << "159";
        //     LoggerMessage message3("Convergence");
        //     message3 << "yuhu";
        //     LoggerMessage message4("Is converged");
        //     message4 << "NO";

        //     p_output->WriteMessage(message1);
        //     p_output->WriteMessage(message2);
        //     p_output->WriteMessage(message3);
        //     p_output->WriteMessage(message4);

        //     // Logger("Time Step") << 123 << std::endl;
        //     // Logger("Iteration Number") << 55 << std::endl;
        //     // Logger("Convergence") << "NO" << std::endl;
        //     // Logger("Is converged") << "Maybe" << std::endl;

        //     std::cout << buffer.str() << std::endl;
        // }
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

        // TODO check if the variables that are being checked here are in the ModelPart!
        // TODO check if the variables that are being checked here Dofs!

        return 0;
        KRATOS_CATCH("");
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GeneralConvergenceCriteria";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GeneralConvergenceCriteria";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    std::vector<TDataType> mRatioTolerances;
    std::vector<TDataType> mAbsTolerances;

    std::vector<TDataType> mRatioResiduals;
    std::vector<TDataType> mAbsResiduals;

    std::vector<SizeType> mNumDofs;

    std::unordered_map<KeyType, IndexType> mKeyToIndexMap;
    std::unordered_map<IndexType, std::string> mIndexToNameMap;

    BasisType mBasisType;


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /*
    Function to print info abt the convergence
    Note that we are flushing the buffer at the end only!
    */
    void PrintConvergenceInfo(const bool IsConverged,
                              const std::vector<bool> ConvergenceInfoVector,
                              const int NonlinIterationNumber)
    {
        const int num_vars_to_separate = mRatioResiduals.size();

        KRATOS_INFO("ConvergenceCriteria") << "Convergence Check; Iteration "
            << NonlinIterationNumber << "\n";

        for (int i=num_vars_to_separate-1; i>-1; --i) // Print in reverse to have the "other dofs" at the end
        {
            std::stringstream conv_info;
            if (ConvergenceInfoVector[i]) conv_info << BOLDFONT(FGRN("converged"));
            else conv_info << BOLDFONT(FRED("not converged"));

            KRATOS_INFO("") << "\t" << mIndexToNameMap[i] << conv_info.str()
                << " | ratio = " << mRatioResiduals[i] << "; exp.ratio = " << mRatioTolerances[i] << " | "
                << "abs = "      << mAbsResiduals[i]   << "; exp.abs = "   << mAbsTolerances[i] << "\n";
        }

        if (IsConverged)
            KRATOS_INFO("ConvergenceCriteria") << BOLDFONT(FGRN("Convergence is achieved in Iteration ")) << NonlinIterationNumber << std::endl; // TODO most likely remove the color
        else
            KRATOS_INFO("ConvergenceCriteria") << BOLDFONT(FRED("Convergence is not achieved in Iteration ")) << NonlinIterationNumber << std::endl; // TODO most likely remove the color
    }


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
    GeneralConvergenceCriteria& operator=(GeneralConvergenceCriteria const& rOther){}

    /// Copy constructor.
    GeneralConvergenceCriteria(GeneralConvergenceCriteria const& rOther){}


    ///@}

}; // Class GeneralConvergenceCriteria

///@}

///@name Type Definitions
///@{


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


