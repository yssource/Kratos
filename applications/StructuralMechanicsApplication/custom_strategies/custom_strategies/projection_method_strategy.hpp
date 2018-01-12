// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Quirin Aumann
//

#if !defined(KRATOS_PROJECTION_METHOD_STRATEGY )
#define  KRATOS_PROJECTION_METHOD_STRATEGY

// System includes
#include<iostream>
#include<vector>
#include<iterator>

// External includes
#include<boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "utilities/qr_utility.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "structural_mechanics_application_variables.h"

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

/// Strategy for solving generalized eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class ProjectionMethodStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ProjectionMethodStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType DenseMatrixPointerType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorType SparseVectorType;
    
    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef std::complex<double> ComplexType;

    typedef boost::numeric::ublas::vector<ComplexType> ComplexVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ProjectionMethodStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver,
        vector< double > SamplingPoints
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpBuilderAndSolver = pBuilderAndSolver;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        mInitializeWasPerformed = false;

        mpForceVector = SparseSpaceType::CreateEmptyVectorPointer();
        mpModalMatrix = DenseSpaceType::CreateEmptyMatrixPointer();

        mpForceVectorReduced = SparseSpaceType::CreateEmptyVectorPointer();
        mpStiffnessMatrixReduced = SparseSpaceType::CreateEmptyMatrixPointer();
        mpMassMatrixReduced = SparseSpaceType::CreateEmptyMatrixPointer();
        mpReducedBasis = SparseSpaceType::CreateEmptyMatrixPointer();

        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;
        mSystemDamping = 0.0;
        mSamplingPoints = SamplingPoints;
        // this->SetUseMaterialDampingFlag(UseMaterialDampingFlag);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build only once)
        this->SetRebuildLevel(0);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    ProjectionMethodStrategy(const ProjectionMethodStrategy& Other) = delete;

    /// Destructor.
    ~ProjectionMethodStrategy() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetIsInitialized(bool val)
    {
        mInitializeWasPerformed = val;
    }

    bool GetIsInitialized() const
    {
        return mInitializeWasPerformed;
    }

    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    };

    SchemePointerType& pGetScheme()
    {
        return mpScheme;
    };

    void SetBuilderAndSolver(BuilderAndSolverPointerType pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    BuilderAndSolverPointerType& pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        this->pGetBuilderAndSolver()->SetReshapeMatrixFlag(flag);
    }

    bool GetReformDofSetAtEachStepFlag() const
    {
        return this->pGetBuilderAndSolver()->GetReshapeMatrixFlag();
    }

    void SetUseMaterialDampingFlag(bool flag)
    {
        mUseMaterialDamping = flag;
    }

    bool GetUseMaterialDampingFlag() const
    {
        return mUseMaterialDamping;
    }

    /// Set verbosity level of the solving strategy.
    /**
     * - 0 -> mute... no echo at all
     * - 1 -> print time and basic information
     * - 2 -> print linear solver data
     * - 3 -> print debug information
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /// Initialization to be performed once before using the strategy.
    virtual void Initialize() override
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        // auto& r_process_info = r_model_part.GetProcessInfo();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Initialize() of ProjectionMethodStrategy." << std::endl;

        this->Check();

        auto& p_scheme = this->pGetScheme();

        if (p_scheme->SchemeIsInitialized() == false)
            p_scheme->Initialize(r_model_part);

        if (p_scheme->ElementsAreInitialized() == false)
            p_scheme->InitializeElements(r_model_part);

        if (p_scheme->ConditionsAreInitialized() == false)
            p_scheme->InitializeConditions(r_model_part);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Initialize() of ProjectionMethodStrategy." << std::endl;

        // set up the system
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();

        // Reset solution dofs
        boost::timer system_construction_time;
        // Set up list of dofs
        boost::timer setup_dofs_time;
        p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;
        }

        // Set global equation ids
        boost::timer setup_system_time;
        p_builder_and_solver->SetUpSystem(r_model_part);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "setup_system_time : " << setup_system_time.elapsed() << std::endl;
        }

        /////////////////////////////////////////////////////
        // get K, M and f
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();

        //initialize dummy vectors
        auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        auto pb = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rDx = *pDx;
        auto& rb = *pb;
        SparseSpaceType::Resize(rDx,system_size);
        SparseSpaceType::Set(rDx,0.0);
        SparseSpaceType::Resize(rb,system_size);
        SparseSpaceType::Set(rb,0.0);

        boost::timer system_build_time;
        
        //mass matrix:
        r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
        
        auto mass_matrix = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_mass_matrix = *mass_matrix;
        p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, 
            mass_matrix,
            pDx,
            pb,
            r_model_part.Elements(),
            r_model_part.Conditions(),
            r_model_part.GetProcessInfo());

        p_builder_and_solver->BuildLHS(p_scheme, r_model_part, r_mass_matrix);
        this->ApplyDirichletConditions(r_mass_matrix, 1.0);

        //stiffness and force
        r_model_part.GetProcessInfo()[BUILD_LEVEL] = 2;
        
        auto stiffness_matrix = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_stiffness_matrix = *stiffness_matrix;
        p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, 
            stiffness_matrix,
            pDx,
            pb,
            r_model_part.Elements(),
            r_model_part.Conditions(),
            r_model_part.GetProcessInfo());

        auto force_vector = SparseSpaceType::CreateEmptyVectorPointer();
        auto& r_force_vector = *force_vector;
        SparseSpaceType::Resize(r_force_vector, system_size);
        SparseSpaceType::Set(r_force_vector,0.0);

        p_builder_and_solver->Build(p_scheme, r_model_part, r_stiffness_matrix, r_force_vector);
        this->ApplyDirichletConditions(r_stiffness_matrix, 1.0);


        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "system_build_time : " << system_build_time.elapsed() << std::endl;
        }

        // KRATOS_WATCH(r_force_vector)
        // KRATOS_WATCH(r_stiffness_matrix)
        // KRATOS_WATCH(r_mass_matrix)

        //sampling points
        KRATOS_WATCH(mSamplingPoints)
        const std::size_t n_sampling_points = mSamplingPoints.size();
        const std::size_t reduced_system_size = 3 * n_sampling_points;

        //initialize sb, As, AAs vectors
        auto s = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rs = *s;
        SparseSpaceType::Resize(rs,system_size);
        SparseSpaceType::Set(rs,0.0);
        auto As = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rAs = *As;
        SparseSpaceType::Resize(rAs,system_size);
        SparseSpaceType::Set(rAs,0.0);
        auto AAs = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rAAs = *AAs;
        SparseSpaceType::Resize(rAAs,system_size);
        SparseSpaceType::Set(rAAs,0.0);

        auto kdyn = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_kdyn = *kdyn;
        SparseSpaceType::Resize(r_kdyn, system_size, system_size);

        auto basis = DenseSpaceType::CreateEmptyMatrixPointer();
        auto r_basis = *basis;
        DenseSpaceType::Resize(r_basis, system_size, reduced_system_size);

        vector< double > aux;
        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            KRATOS_WATCH( mSamplingPoints(i) )
            // KRATOS_WATCH( std::pow( mSamplingPoints(i), 2.0) )
            r_kdyn = r_stiffness_matrix - ( std::pow( mSamplingPoints(i), 2.0 ) * r_mass_matrix );
            // KRATOS_WATCH(r_kdyn)
            // KRATOS_WATCH(rs)
            // KRATOS_WATCH(r_force_vector)
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, rs, r_force_vector );
            // KRATOS_WATCH(rs)
            aux = prod( r_mass_matrix, rs );
            // KRATOS_WATCH(aux)
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, rAs, aux );
            aux = prod( r_mass_matrix, rAs );
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, rAAs, aux );

            // KRATOS_WATCH(rs)
            // KRATOS_WATCH(rAs)
            // KRATOS_WATCH(rAAs)

            column( r_basis, (i*3) ) = rs;
            column( r_basis, (i*3)+1 ) = rAs;
            column( r_basis, (i*3)+2 ) = rAAs;


            // KRATOS_WATCH(r_basis)


        }
        // DenseMatrixPointerType pAuxMatQR(new DenseMatrixType(system_size, 3*n_sampling_points));
        // auto& r_pAuxMatQR = *pAuxMatQR;
        // auto ff = (r_pAuxMatQR)(0,0);
        // auto* la = &(r_pAuxMatQR)(0,0);
        // std::cout << "yo" << std::endl;
        // KRATOS_WATCH(&(r_basis)(0,0))

        //orthogonalize the basis -> basis_r
        mQR_decomposition.compute( system_size, 3*n_sampling_points, &(r_basis)(0,0) );
        // std::cout << "yo2" << std::endl;
        // KRATOS_WATCH(r_basis)
        mQR_decomposition.compute_q();

        // auto basis_r = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_basis_r = *mpReducedBasis;
        SparseSpaceType::Resize(r_basis_r, system_size, 3*n_sampling_points);
        for( size_t i = 0; i < system_size; ++i )
        {
            for( size_t j = 0; j < (3*n_sampling_points); ++j )
            {
                r_basis_r(i,j) = mQR_decomposition.Q(i,j);
            }
        }
        // KRATOS_WATCH(r_basis)
        // KRATOS_WATCH(r_basis_r)

        // auto mat_r = SparseSpaceType::CreateEmptyMatrixPointer();
        // auto& r_mat_r = *mat_r;
        // SparseSpaceType::Resize(r_mat_r, 3*n_sampling_points, 3*n_sampling_points);
        // for( size_t i = 0; i < 3*n_sampling_points; ++i )
        // {
        //     for( size_t j = 0; j < (3*n_sampling_points); ++j )
        //     {
        //         r_mat_r(i,j) = mQR_decomposition.R(i,j);
        //     }
        // }
        // KRATOS_WATCH(r_mat_r)

        //reduce mass, stiffness, force
        // mpMassMatrixReduced = prod( prod( r_basis_r, mass_matrix ), r_basis_r );
        // mpStiffnessMatrixReduced = prod( prod( r_basis_r, stiffness_matrix ), r_basis_r );
        // auto hhh = prod( r_basis_r, r_force_vector );
        // KRATOS_WATCH(prod( r_force_vector, r_basis_r ))
        // mpForceVectorReduced = hhh;
        auto& r_force_vector_reduced = *mpForceVectorReduced;
        // mpForceVectorReduced = ZeroMatrix( 3*n_sampling_points );
        r_force_vector_reduced = (prod( r_force_vector, r_basis_r ));
        auto& r_stiffness_matrix_reduced = *mpStiffnessMatrixReduced;
        auto& r_mass_matrix_reduced = *mpMassMatrixReduced;

        // KRATOS_WATCH( prod( matrix<double>(prod(trans(r_basis_r),r_stiffness_matrix)), r_basis_r))
        r_stiffness_matrix_reduced = prod( matrix< double >( prod( trans( r_basis_r ),r_stiffness_matrix ) ), r_basis_r );
        r_mass_matrix_reduced = prod( matrix< double >( prod( trans( r_basis_r ),r_mass_matrix ) ), r_basis_r );
        
        // KRATOS_WATCH( prod( aa,r_basis_r ) )
        // KRATOS_WATCH( prod( ))

        // KRATOS_WATCH( prod( prod( trans(r_basis_r), r_stiffness_matrix ), r_basis_r ) )
        

        // KRATOS_WATCH(*mpMassMatrixReduced)
        // KRATOS_WATCH(*mpStiffnessMatrixReduced)
        // KRATOS_WATCH(*mpForceVectorReduced)
        std::cout << "ende" << std::endl;
        
        
        KRATOS_CATCH("")
    }

    double Solve() override
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        // std::cout << "****************************++solve***********************+" << std::endl;
        // operations to be done once
        if (this->GetIsInitialized() == false)
        {
            Initialize();
            this->SetIsInitialized(true);
        }

        this->InitializeSolutionStep();

        auto& r_process_info = r_model_part.GetProcessInfo();
        double excitation_frequency = r_process_info[TIME];
        const unsigned int system_size = this->pGetBuilderAndSolver()->GetEquationSystemSize();
        const std::size_t reduced_system_size = 3 * mSamplingPoints.size();

        auto& r_stiffness_matrix_reduced = *mpStiffnessMatrixReduced;
        auto& r_mass_matrix_reduced = *mpMassMatrixReduced;
        auto& r_force_vector_reduced = *mpForceVectorReduced;

        // auto displacement_reduced = SparseSpaceType::CreateEmptyVectorPointer();
        // auto& r_displacement_reduced = *displacement_reduced;
        // SparseSpaceType::Resize(r_displacement_reduced,reduced_system_size);
        // SparseSpaceType::Set(r_displacement_reduced,0.0);

        SparseVectorType displacement_reduced;
        displacement_reduced.resize( reduced_system_size, false );
        displacement_reduced = ZeroVector( reduced_system_size );

        // KRATOS_WATCH(excitation_frequency)
        // KRATOS_WATCH(r_stiffness_matrix_reduced)
        // KRATOS_WATCH(r_mass_matrix_reduced)

        auto kdyn = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_kdyn = *kdyn;
        SparseSpaceType::Resize(r_kdyn, reduced_system_size, reduced_system_size);

        r_kdyn = r_stiffness_matrix_reduced - ( std::pow( excitation_frequency, 2.0 ) * r_mass_matrix_reduced );
        // KRATOS_WATCH(r_kdyn)

        this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, displacement_reduced, r_force_vector_reduced );

        SparseVectorType displacement;
        displacement.resize( system_size, false );
        displacement = ZeroVector( system_size );

        displacement = prod( *mpReducedBasis, displacement_reduced );
        // KRATOS_WATCH(displacement_reduced)
        // KRATOS_WATCH(displacement)
        

        this->AssignVariables(displacement);
        this->FinalizeSolutionStep();

        return 0.0;

        KRATOS_CATCH("")
    }

    /// Clear the strategy.
    virtual void Clear() override
    {
        KRATOS_TRY

        // if the preconditioner is saved between solves, it should be cleared here
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();
        p_builder_and_solver->GetLinearSystemSolver()->Clear();

        SparseSpaceType::Clear(mpForceVector);
        DenseSpaceType::Clear(mpModalMatrix);

        // re-setting internal flag to ensure that the dof sets are recalculated
        p_builder_and_solver->SetDofSetIsInitializedFlag(false);

        p_builder_and_solver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;
        mUseMaterialDamping = false;
        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;
        mSystemDamping = 0.0;

        KRATOS_CATCH("")
    }

    /// Initialization to be performed before every solve.
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering InitializeSolutionStep() of ProjectionMethodStrategy" << std::endl;

        BuilderAndSolverPointerType& p_builder_and_solver = this->pGetBuilderAndSolver();
        SchemePointerType& p_scheme = this->pGetScheme();
        auto& r_force_vector = *mpForceVector;

        // // initialize dummy vectors
        auto pA = SparseSpaceType::CreateEmptyMatrixPointer();
        auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rA = *pA;
        auto& rDx = *pDx;

        SparseSpaceType::Resize(rA,SparseSpaceType::Size(r_force_vector),SparseSpaceType::Size(r_force_vector));
        SparseSpaceType::SetToZero(rA);
        SparseSpaceType::Resize(rDx,SparseSpaceType::Size(r_force_vector));
        SparseSpaceType::Set(rDx,0.0);

        // initial operations ... things that are constant over the solution step
        p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,r_force_vector);

        // initial operations ... things that are constant over the solution step
        p_scheme->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,r_force_vector);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting InitializeSolutionStep() of ProjectionMethodStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /// Check whether initial input is valid.
    virtual int Check() override
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Check() of ProjectionMethodStrategy" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(r_model_part);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(r_model_part);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Check() of ProjectionMethodStrategy" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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

    ///@}
    ///@name Member Variables
    ///@{

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    bool mInitializeWasPerformed;

    SparseVectorPointerType mpForceVector;

    SparseVectorPointerType mpForceVectorReduced;

    SparseMatrixPointerType mpMassMatrixReduced;

    SparseMatrixPointerType mpStiffnessMatrixReduced;

    SparseMatrixPointerType mpReducedBasis;

    DenseMatrixPointerType mpModalMatrix;

    double mRayleighAlpha;

    double mRayleighBeta;

    double mSystemDamping;

    bool mUseMaterialDamping;

    vector< double > mMaterialDampingRatios;

    vector< double > mSamplingPoints;

    QR<double, row_major>                 mQR_decomposition;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Assign the modal displacement to the dofs and the phase angle to the reaction
    void AssignVariables(SparseVectorType& rDisplacement, int step=0)
    {
        auto& r_model_part = BaseType::GetModelPart();
        for( auto& node : r_model_part.Nodes() )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = node.GetDofs();
            
            for( auto it_dof = std::begin(rNodeDofs); it_dof != std::end(rNodeDofs); it_dof++ )
            {
                if( !it_dof->IsFixed() )
                {
                    //absolute displacement
                    it_dof->GetSolutionStepValue(step) = rDisplacement(it_dof->EquationId());
                    //phase angle
                    // it_dof->GetSolutionStepReactionValue(step) = std::abs(std::arg(rModalDisplacement(it_dof->EquationId())));
                }
                else
                {
                    it_dof->GetSolutionStepValue(step) = 0.0;
                }
            }
        }
    }

    /// Apply Dirichlet boundary conditions without modifying dof pattern.
    /**
     *  The dof pattern is preserved to support algebraic multigrid solvers with
     *  component-wise aggregation. Rows and columns of the fixed dofs are replaced
     *  with zeros on the off-diagonal and the diagonal is scaled by factor.
     *  Taken from eigensolver_strategy.hpp
     */
    void ApplyDirichletConditions(
        SparseMatrixType& rA, 
        double Factor
        )
    {
        KRATOS_TRY

        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
        {
            std::cout << "Entering ApplyDirichletConditions() of EigensolverStrategy" << std::endl;
        }

        const std::size_t SystemSize = rA.size1();
        std::vector<double> ScalingFactors(SystemSize);
        auto& rDofSet = this->pGetBuilderAndSolver()->GetDofSet();
        const int NumDofs = static_cast<int>(rDofSet.size());

        // NOTE: dofs are assumed to be numbered consecutively
        #pragma omp parallel for firstprivate(NumDofs)
        for(int k = 0; k<NumDofs; k++)
        {
            auto dof_iterator = std::begin(rDofSet) + k;
            ScalingFactors[k] = (dof_iterator->IsFixed()) ? 0.0 : 1.0;
        }

        double* AValues = std::begin(rA.value_data());
        std::size_t* ARowIndices = std::begin(rA.index1_data());
        std::size_t* AColIndices = std::begin(rA.index2_data());

        // if there is a line of all zeros, put one on the diagonal
        // #pragma omp parallel for firstprivate(SystemSize)
        // for(int k = 0; k < static_cast<int>(SystemSize); ++k)
        // {
        //     std::size_t ColBegin = ARowIndices[k];
        //     std::size_t ColEnd = ARowIndices[k+1];
        //     bool empty = true;
        //     for (auto j = ColBegin; j < ColEnd; ++j)
        //         if(AValues[j] != 0.0)
        //         {
        //             empty = false;
        //             break;
        //         }
        //     if(empty == true)
        //         rA(k,k) = 1.0;
        // }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(SystemSize); ++k)
        {
            std::size_t ColBegin = ARowIndices[k];
            std::size_t ColEnd = ARowIndices[k+1];
            if (ScalingFactors[k] == 0.0)
            {
                // row dof is fixed. zero off-diagonal columns and factor diagonal
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                {
                    if (static_cast<int>(AColIndices[j]) != k)
                    {
                        AValues[j] = 0.0;
                    }
                    else
                    {
                        AValues[j] *= Factor;
                    }
                }
            }
            else
            {
                // row dof is not fixed. zero columns associated with fixed dofs
                for (std::size_t j = ColBegin; j < ColEnd; ++j)
                {
                    AValues[j] *= ScalingFactors[AColIndices[j]];
                }
            }
        }

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
        {
            std::cout << "Exiting ApplyDirichletConditions() of EigensolverStrategy" << std::endl;
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class ProjectionMethodStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos */

#endif /* KRATOS_PROJECTION_METHOD_STRATEGY  defined */

