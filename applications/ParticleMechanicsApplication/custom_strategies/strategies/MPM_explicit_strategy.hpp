//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_EXPLICIT_STRATEGY )
#define  KRATOS_MPM_EXPLICIT_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "includes/kratos_flags.h"

#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

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

/// Short class definition.

/**   Detail class definition.



 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class MPMExplicitStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(MPMExplicitStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::ElementsArrayType ElementsArrayType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructors.
     */
    MPMExplicitStrategy(
        ModelPart& model_part,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
    }


    MPMExplicitStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver)
                             );

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;

        mInitializeWasPerformed = false;

        mFinalizeSolutionStep = true;

        // Tells to the Builder And Solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH( "" )
    }

    MPMExplicitStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY

        mKeepSystemConstantDuringIterations = false;

        // Set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;

        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        // Saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        // Saving the scheme
        mpScheme = pScheme;

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        // Setting up the default builder and solver
        mpBuilderAndSolver = pNewBuilderAndSolver;

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;

        mInitializeWasPerformed = false;

        mFinalizeSolutionStep = true;

        // Tells to the Builder And Solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH( "" )
    }

    /** Destructor.
     */
    virtual ~MPMExplicitStrategy()
    {
    }

    /** Destructor.
     */

    //Set and Get Scheme ... containing Builder, Update and other

    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    //Set and Get the BuilderAndSolver

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        mReformDofSetAtEachStep = flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    void SetFinalizeSolutionStepFlag(bool FinalizeSolutionStepFlag = true)
    {
        mFinalizeSolutionStep = FinalizeSolutionStepFlag;
    }

    bool GetFinalizeSolutionStepFlag()
    {
        return mFinalizeSolutionStep;
    }

    //level of echo for the solving strategy
    // 0 -> mute... no echo at all
    // 1 -> printing time and basic informations
    // 2 -> printing linear solver data
    // 3 -> Print of debug informations:
    //		Echo of stiffness matrix, Dx, b...

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        /* GetBuilderAndSolver()->SetEchoLevel(Level);
        mpConvergenceCriteria->SetEchoLevel(Level); */
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
    //*********************************************************************************
    /**
    Initialize members
     */
    //**********************************************************************


    void Initialize() override
    {
        KRATOS_TRY

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();

        // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        // if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
        {
            KRATOS_INFO_IF("MPM_Strategy",this->GetEchoLevel() >1) << "Initializing solving strategy" << std::endl;
            KRATOS_ERROR_IF(mInitializeWasPerformed == true) << "Initialization was already performed " << mInitializeWasPerformed << std::endl;

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            KRATOS_INFO_IF("MPM_Strategy",this->GetEchoLevel() >1) << "Initializing scheme" << std::endl;
            if (pScheme->SchemeIsInitialized() == false)
                pScheme->Initialize(BaseType::GetModelPart());

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            KRATOS_INFO_IF("MPM_Strategy",this->GetEchoLevel() >1) << "Initializing elements" << std::endl;
            if (pScheme->ElementsAreInitialized() == false)
                pScheme->InitializeElements(BaseType::GetModelPart());

            // Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            KRATOS_INFO_IF("MPM_Strategy",this->GetEchoLevel() >1) << "Initializing conditions" << std::endl;
            if (pScheme->ConditionsAreInitialized() == false)
                pScheme->InitializeConditions(BaseType::GetModelPart());


            /* //todo: delete this
            // Initialisation of the convergence criteria
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
            KRATOS_INFO_IF("MPM_Strategy",this->GetEchoLevel() >1) << "Initializing convergence criteria"<<std::endl;
            if (mpConvergenceCriteria->IsInitialized() == false)
                mpConvergenceCriteria->Initialize(BaseType::GetModelPart()); */

            mInitializeWasPerformed = true;
        }

        /* // Set up the system, operation performed just once unless it is required to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
        {
            // Setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme, BaseType::GetModelPart());

            // Shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());

        } */

        // Prints informations about the current time
        if (this->GetEchoLevel() == 2 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            KRATOS_INFO("MPM_Strategy") << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
        }

        KRATOS_CATCH( "" )
    }


    //*********************************************************************************
    /**
    the problem of interest is solved
     */
    //**********************************************************************
    bool SolveSolutionStep() override
    {
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        DofsArrayType dof_set_dummy;
        TSystemMatrixType mA = TSystemMatrixType();
        TSystemVectorType mDx = TSystemVectorType();
        TSystemVectorType mb = TSystemVectorType();


        // Initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        // Compute residual forces on the model part
        this->CalculateAndAddRHS(pScheme, BaseType::GetModelPart());

        pScheme->Update(BaseType::GetModelPart(), dof_set_dummy, mA, mDx, mb);

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        // Calculate reactions if required
        if (mCalculateReactionsFlag) {
            CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }

        return true;
    }

    //*********************************************************************************

    //**********************************************************************
    //**********************************************************************

    void Clear() override
    {
        KRATOS_TRY
        /* // Setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear(); */

        GetScheme()->Clear();

        KRATOS_CATCH( "" )
    }

    /*@} */
    /**@name Operators
     */
    /*@{ */

    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */

    /*@{ */

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }

    void SetKeepSystemConstantDuringIterations(bool value)
    {
        mKeepSystemConstantDuringIterations = value;
    }

    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
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

protected:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    typename TSchemeType::Pointer mpScheme;

    typename TLinearSolver::Pointer mpLinearSolver;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;

    /*		TSystemVectorType mDx;
                    TSystemVectorType mb;
                    TSystemMatrixType mA;*/
    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;

    /**
    Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    - true  => reforme at each time step
    - false => form just one (more efficient)

    Default = false
     */
    bool mReformDofSetAtEachStep;

    /**
    Flag telling if it is needed or not to compute the reactions

    default = true
     */
    bool mCalculateReactionsFlag;

    bool mSolutionStepIsInitialized;

    ///default = 30
    unsigned int mMaxIterationNumber;

    bool mInitializeWasPerformed;

    //flag to allow keeping system matrix constant during iterations
    bool mKeepSystemConstantDuringIterations;

    //flag to allow to not finalize the solution step, so the historical variables are not updated
    bool mFinalizeSolutionStep;

    /*@} */
    /**@name Private Operators*/
    /*@{ */
    //**********************************************************************
    //**********************************************************************

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        // Initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
            typename TSchemeType::Pointer pScheme = GetScheme();

            ModelPart& r_model_part = BaseType::GetModelPart();

            TSystemMatrixType matrix_a_dummy = TSystemMatrixType();
            TSystemVectorType rDx = TSystemVectorType();
            TSystemVectorType rb = TSystemVectorType();

            // Setting up the Vectors involved to the correct size
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA, mpDx, mpb, BaseType::GetModelPart());

            TSystemMatrixType mA = TSystemMatrixType();
            TSystemVectorType mDx = TSystemVectorType();
            TSystemVectorType mb = TSystemVectorType();

            // Initial operations ... things that are constant over the Solution Step
            pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

            // Initial operations ... things that are constant over the Solution Step
            pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);


            if (BaseType::mRebuildLevel > 0)
            { // TODO: Right now is computed in the Initialize() because is always zero, the option to set the RebuildLevel should be added in the constructor or in some place
                ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
                ElementsArrayType& r_elements = r_model_part.Elements();
                const auto it_elem_begin = r_elements.begin();

                // Set Nodal Mass and Damping to zero
                NodesArrayType& r_nodes = r_model_part.Nodes();
                VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
                VariableUtils().SetNonHistoricalVariable(NODAL_DISPLACEMENT_DAMPING, 0.0, r_nodes);

                Vector dummy_vector;
                // If we consider the rotation DoF
                const bool has_dof_for_rot_z = r_model_part.Nodes().begin()->HasDofFor(ROTATION_Z);
                if (has_dof_for_rot_z) {
                    const array_1d<double, 3> zero_array = ZeroVector(3);
                    VariableUtils().SetNonHistoricalVariable(NODAL_INERTIA, zero_array, r_nodes);
                    VariableUtils().SetNonHistoricalVariable(NODAL_ROTATION_DAMPING, zero_array, r_nodes);

                    #pragma omp parallel for firstprivate(dummy_vector), schedule(guided,512)
                    for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                        // Getting nodal mass and inertia from element
                        // this function needs to be implemented in the respective
                        // element to provide inertias and nodal masses
                        auto it_elem = it_elem_begin + i;
                        it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_INERTIA, r_current_process_info);
                    }
                } else { // Only NODAL_MASS and NODAL_DISPLACEMENT_DAMPING are needed
                    #pragma omp parallel for firstprivate(dummy_vector), schedule(guided,512)
                    for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                        // Getting nodal mass and inertia from element
                        // this function needs to be implemented in the respective
                        // element to provide nodal masses
                        auto it_elem = it_elem_begin + i;
                        it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
                    }
            }
            }
            mSolutionStepIsInitialized = true;
        }

        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() >= 3) << "Initialize Solution Step in strategy finished" <<std::endl;

        KRATOS_CATCH( "" )
    }


    //**********************************************************************
    //**********************************************************************
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = GetScheme();
        TSystemMatrixType mA = TSystemMatrixType();
        TSystemVectorType mDx = TSystemVectorType();
        TSystemVectorType mb = TSystemVectorType();
        /* if (mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        } */

        // Calling rDofSet
        //DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        /*Finalization of the solution step,
        operations to be done after achieving convergence, for example the
        Final Residual Vector (mb) has to be saved in there
        to avoid error accumulation*/
        if( mFinalizeSolutionStep )
        {
            KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() >= 3) << "Calling FinalizeSolutionStep" <<std::endl;

            pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
            //pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
            //mpConvergenceCriteria->FinalizeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
        }

        // Cleaning memory after the solution
        pScheme->Clean();

        // Reset flags for next step
        mSolutionStepIsInitialized = false;
        KRATOS_CATCH( "" )
    }

    /**
     * function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();
        GetScheme()->Check(BaseType::GetModelPart());
        return 0;

        KRATOS_CATCH( "" )
    }


    void CalculateAndAddRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        ConditionsArrayType& r_conditions = rModelPart.Conditions();
        ElementsArrayType& r_elements = rModelPart.Elements();

        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
        Element::EquationIdVectorType equation_id_vector_dummy; // Dummy

        #pragma omp parallel for firstprivate(RHS_Contribution, equation_id_vector_dummy), schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_conditions.size()); ++i) {
            auto it_cond = r_conditions.begin() + i;
            pScheme->Condition_Calculate_RHS_Contribution((*it_cond.base()), RHS_Contribution, equation_id_vector_dummy, r_current_process_info);
        }

        #pragma omp parallel for firstprivate(RHS_Contribution, equation_id_vector_dummy), schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
            auto it_elem = r_elements.begin() + i;
            pScheme->Calculate_RHS_Contribution((*it_elem.base()), RHS_Contribution, equation_id_vector_dummy, r_current_process_info);
        }

        KRATOS_CATCH("")
    }


    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        // We iterate over the nodes
        auto& r_nodes = rModelPart.Nodes();

        // If we consider rotation dofs
        const bool has_dof_for_rot_z = (r_nodes.begin())->HasDofFor(ROTATION_Z);

        // Auxiliar values
        const array_1d<double, 3> zero_array = ZeroVector(3);
        array_1d<double, 3> force_residual = ZeroVector(3);
        array_1d<double, 3> moment_residual = ZeroVector(3);

        // Getting
        const auto it_node_begin = r_nodes.begin();
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
        const IndexType rotppos = it_node_begin->GetDofPosition(ROTATION_X);

        // Iterating nodes
        #pragma omp parallel for firstprivate(force_residual, moment_residual), schedule(guided,512)
        for(int i=0; i<static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = it_node_begin + i;

            noalias(force_residual) = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
            if (has_dof_for_rot_z) {
                noalias(moment_residual) = it_node->FastGetSolutionStepValue(MOMENT_RESIDUAL);
            } else {
                noalias(moment_residual) = zero_array;
            }

            if (it_node->GetDof(DISPLACEMENT_X, disppos).IsFixed()) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_X);
                r_reaction = force_residual[0];
            }
            if (it_node->GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed()) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Y);
                r_reaction = force_residual[1];
            }
            if (it_node->GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed()) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Z);
                r_reaction = force_residual[2];
            }
            if (has_dof_for_rot_z) {
                if (it_node->GetDof(ROTATION_X, rotppos).IsFixed()) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_MOMENT_X);
                    r_reaction = moment_residual[0];
                }
                if (it_node->GetDof(ROTATION_Y, rotppos + 1).IsFixed()) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_MOMENT_Y);
                    r_reaction = moment_residual[1];
                }
                if (it_node->GetDof(ROTATION_Z, rotppos + 2).IsFixed()) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_MOMENT_Z);
                    r_reaction = moment_residual[2];
                }
            }
        }
    }



    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /** Copy constructor.
     */
    MPMExplicitStrategy(const MPMExplicitStrategy& Other)
    {
    };


    /*@} */

}; /* Class MPMExplicitStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}; /* namespace Kratos.*/

#endif /* KRATOS_MPM_EXPLICIT_STRATEGY  defined */

