// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Aditya Ghantasala
//

#if !defined(LINEAR_MOR_MATRIX_OUTPUT_STRATEGY )
#define  LINEAR_MOR_MATRIX_OUTPUT_STRATEGY

/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>

/* External includes */

/* Project includes */
// #include "structural_mechanics_application.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h" // TODO: change to new b&s

#include <cmath>

namespace Kratos
{
template<class TSparseSpace,
         class TDenseSpace, 
         class TLinearSolver 
         >
class LinearMorMatrixOutputStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

    typedef LineSearchesUtility<TSparseSpace, TDenseSpace, TLinearSolver> TlineSearchesType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( LinearMorMatrixOutputStrategy );

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

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

    typedef long double RealType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /************************************* CONSTRUCTOR *********************************/
    /***********************************************************************************/
    
    LinearMorMatrixOutputStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            bool MoveMeshFlag           = true
            )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
    {
        KRATOS_TRY;

        // Saving the scheme
        mpScheme = pScheme;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
                             );

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed    = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // Be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(false);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH("");
    }

    /************************************* DESTRUCTOR **********************************/
    /***********************************************************************************/
    
    ~ResidualBasedArcLengthStrategy() override {}

    /************************************* OPERATIONS **********************************/
    /***********************************************************************************/

    //Set and Get Scheme ... containing Builder, Update and other
    void SetScheme(typename TSchemeType::Pointer pScheme )
    {
        mpScheme = pScheme;
    };
    
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    // Set and Get the BuilderAndSolver
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver )
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };
    
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

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

    void SetMaxIterationNumber(unsigned int  MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }
    
    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    // Level of echo for the solving strategy
    // 0 -> Mute... no echo at all
    // 1 -> Printing time and basic informations
    // 2 -> Printing linear solver data
    // 3 -> Print of debug informations:
    // Echo of stiffness matrix, Dx, b...
    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    * values of the solution step of interest are assumed equal to the old values
    */
    
    void Predict() override
    {
        KRATOS_TRY;
        
        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& mA  = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb  = *mpb;
	
        GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);

        KRATOS_CATCH("");
    }


    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It solves the problem
    */
    
    double Solve() override
    {
        KRATOS_TRY;

        //std::cout<<std::fixed<<std::setw(15)<<std::scientific<<std::setprecision(9);
        if (this->GetEchoLevel() > 0)
        {
            std::cout<<"************************************************************************"<<std::endl;
            std::cout<<"Begininning Linear MOR strategy to print out matrices ... "<<std::endl;
            std::cout<<"************************************************************************"<<std::endl;
        }


        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        //ModelPart& r_model_part = BaseType::GetModelPart();
	
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
	
        // Creating models part for analysis
        InitializeAuxiliaryModelParts(BaseType::GetModelPart());
        mstep = BaseType::GetModelPart().GetProcessInfo()[STEP];

        // Initialisation of the convergence criteria and variables of arc lenght
        if(mInitializeWasPerformed == false)
        {
            Initialize();
        }

        // Set up the system, operation performed just once unless it is required to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true )
        {
            // Setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme,BaseType::GetModelPart());

            // Shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
        }
        
        // Updates the database with a prediction of the solution
        Predict();

        // Initialize solution step
        if (mSolutionStepIsInitialized == false)
        {
            InitializeSolutionStep();
        }


        pBuilderAndSolver->BuildRHS(pScheme,mAuxElementModelPart, mb);

        TSparseSpace::SetToZero(mDx);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mA);
        pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);

        pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);

    }
	

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It clears the variables of the arc length
    */
    
    void Clear() override
    {
        KRATOS_TRY;
        if (this->GetEchoLevel() > 0)
        {
            std::cout << "Arc Length Strategy  Clear function used" << std::endl;
        }

        TSystemMatrixType& mA          = *mpA;
        TSystemVectorType& mDx         = *mpDx;
        TSystemVectorType& mb          = *mpb;
        TSystemVectorType& mDelta_p    = *mpDelta_p;
        TSystemVectorType& mDelta_pold = *mpDelta_pold;

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Resize(mA, 0, 0);

        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Resize(mDx, 0);

        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Resize(mb, 0);

        SparseSpaceType::Clear(mpDelta_p);
        SparseSpaceType::Resize(mDelta_p, 0);

        SparseSpaceType::Clear(mpDelta_pold);
        SparseSpaceType::Resize(mDelta_pold, 0);
	
        // Setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();

        KRATOS_CATCH("");
    }

    /**
    * Returns the LHS of the problem
    */

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }

    /***********************************************************************************/
    /***********************************************************************************/

protected:

private:

    typename TSchemeType::Pointer mpScheme;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    TSystemVectorPointerType mpRHS;
    TSystemMatrixPointerType mpA;
    TSystemMatrixPointerType mpM;


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
    bool mInitializeWasPerformed;

    bool mSolutionStepIsInitialized;
    unsigned int mMaxIterationNumber;
    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Initilise the variables, schemes and convergence criterias
    */

    void Initialize() override
    {
        KRATOS_TRY;

        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
        if (pScheme->SchemeIsInitialized() == false)
        {
            pScheme->Initialize(BaseType::GetModelPart());
        }

        // Initialize The Elements - OPERATIONS TO BE DONE ONCE
        if (pScheme->ElementsAreInitialized() == false)
        {
            pScheme->InitializeElements(BaseType::GetModelPart());
        }

        // Initialize Conditions
        if (pScheme->ConditionsAreInitialized() == false)
        {
            pScheme->InitializeConditions(mAuxConditionModelPart);
        }

        // Initialisation of the convergence criteria
        if (mpConvergenceCriteria->IsInitialized() == false)
        {
            pConvergenceCriteria->Initialize(BaseType::GetModelPart());
        }

        mInitializeWasPerformed = true;

        VariablesArcLength(); // Initialising the variables of the arc length

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It initialises the solution step
    */

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->GetEchoLevel() > 0)
        {
            std::cout<< "Initializing Solution Step " << std::endl;
        }
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = GetScheme();

        // Setting up the Vectors involved to the correct size with value cero
        pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA,mpDx,mpb,BaseType::GetModelPart());

        TSystemMatrixType& mA            = *mpA;          // Stiffness Matrix
        TSystemMatrixType& mM            = *mpM;           // Mass Matrix
        TSystemVectorType& mRHS          = *mpRHS;          // RHS vector
	
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mM);
        TSparseSpace::SetToZero(mRHS);

        // Initial operations
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mRHS,mRHS);
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mRHS,mRHS);

        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mM,mRHS,mRHS);
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mM,mRHS,mRHS);

        mSolutionStepIsInitialized = true;

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * It finalises the arc length for the currrent step
    * @param mIterationNumber: The iteration number in the non-linear step
    * @param mReduceArcLenght: Boolean that tells if the arc length has been computed with the reduced method
    */

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver       = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme                           = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
	
        TSystemMatrixType& mA            = *mpA;
        TSystemMatrixType& mM            = *mpM;
        TSystemVectorType& mRHS            = *mpRHS;

        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mRHS,mRHS);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mRHS,mRHS);

        pScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mM,mRHS,mRHS);
        pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mM,mRHS,mRHS);

        KRATOS_CATCH("");
    }

    LinearMorMatrixOutputStrategy(const LinearMorMatrixOutputStrategy& Other);

    /*@} */

}; /* Class LinearMorMatrixOutputStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* MOR_MATRIX_OUTPUT_STRATEGY  defined */

