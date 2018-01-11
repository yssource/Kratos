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

        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "system_build_time : " << system_build_time.elapsed() << std::endl;
        }

        KRATOS_WATCH(r_force_vector)
        KRATOS_WATCH(r_stiffness_matrix)
        KRATOS_WATCH(r_mass_matrix)

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
            // KRATOS_WATCH( mSamplingPoints(i) )
            // KRATOS_WATCH( std::pow( mSamplingPoints(i), 2.0) )
            r_kdyn = r_stiffness_matrix - ( std::pow( mSamplingPoints(i), 2.0 ) * r_mass_matrix );
            KRATOS_WATCH(r_kdyn)
            KRATOS_WATCH(rs)
            KRATOS_WATCH(r_force_vector)
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, rs, r_force_vector );
            KRATOS_WATCH(rs)
            aux = prod( r_mass_matrix, rs );
            // KRATOS_WATCH(aux)
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, rAs, aux );
            aux = prod( r_mass_matrix, rAs );
            this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve( r_kdyn, rAAs, aux );

            KRATOS_WATCH(rs)
            KRATOS_WATCH(rAs)
            KRATOS_WATCH(rAAs)

            column( r_basis, (i*3) ) = rs;
            column( r_basis, (i*3)+1 ) = rAs;
            column( r_basis, (i*3)+2 ) = rAAs;


            KRATOS_WATCH(r_basis)


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
        KRATOS_WATCH(r_basis)
        mQR_decomposition.compute_q();

        auto basis_r = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_basis_r = *basis_r;
        SparseSpaceType::Resize(r_basis_r, system_size, 3*n_sampling_points);
        for( size_t i = 0; i < system_size; ++i )
        {
            for( size_t j = 0; j < (3*n_sampling_points); ++j )
            {
                r_basis_r(i,j) = mQR_decomposition.Q(i,j);
            }
        }
        KRATOS_WATCH(r_basis_r)

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
        KRATOS_WATCH(prod( r_force_vector, r_basis_r ))
        // mpForceVectorReduced = hhh;
        auto& r_force_vector_reduced = *mpForceVectorReduced;
        // mpForceVectorReduced = ZeroMatrix( 3*n_sampling_points );
        r_force_vector_reduced = (prod( r_force_vector, r_basis_r ));
        auto& r_stiffness_matrix_reduced = *mpStiffnessMatrixReduced;
        // KRATOS_WATCH(prod(r_stiffness_matrix, r_basis_r))
        // auto aa = prod(r_stiffness_matrix, r_basis_r);
        // auto aaa = prod( r_basis_r, trans( aa ) );
        // KRATOS_WATCH( prod( prod( r_basis_r, r_mass_matrix ), r_basis_r ) )
        // r_stiffness_matrix_reduced = prod( prod( r_stiffness_matrix, r_basis_r ), r_basis_r );
        auto& r_mass_matrix_reduced = *mpMassMatrixReduced;
        // r_mass_matrix_reduced = prod( prod( r_mass_matrix, r_basis_r ), r_basis_r );
        


        // KRATOS_WATCH( prod( r_basis_r, r_stiffness_matrix ) )
        // KRATOS_WATCH( prod( trans(r_basis_r), r_stiffness_matrix ) )
        // KRATOS_WATCH( prod( r_stiffness_matrix, r_basis_r ) )
        // KRATOS_WATCH( prod( r_stiffness_matrix, trans(r_basis_r) ) )

        // auto aa = prod( trans(r_basis_r), r_stiffness_matrix, r_basis_r );
        std::cout << "Ü************************" << std::endl;
        // KRATOS_WATCH(aa)
        KRATOS_WATCH( prod( matrix<double>(prod(trans(r_basis_r),r_stiffness_matrix)), r_basis_r))
        r_stiffness_matrix_reduced = prod( matrix< double >( prod( trans( r_basis_r ),r_stiffness_matrix ) ), r_basis_r );
        r_mass_matrix_reduced = prod( matrix< double >( prod( trans( r_basis_r ),r_mass_matrix ) ), r_basis_r );
        
        // KRATOS_WATCH( prod( aa,r_basis_r ) )
        // KRATOS_WATCH( prod( ))

        // KRATOS_WATCH( prod( prod( trans(r_basis_r), r_stiffness_matrix ), r_basis_r ) )
        

        KRATOS_WATCH(*mpMassMatrixReduced)
        KRATOS_WATCH(*mpStiffnessMatrixReduced)
        KRATOS_WATCH(*mpForceVectorReduced)

        
        // mQR_decomposition.compute(system_size, 3*n_sampling_points, &(*pAuxMatQR)(0,0));

            // this->pGetBuilderAndSolver()->GetLinearSystemSolver()->Solve(
            //     rStiffnessMatrix,
            //     rMassMatrix,
            //     Eigenvalues,
            //     Eigenvectors);
        /////////////////////////////////////////////////////

        // initialize the force vector; this does not change during the computation
        // auto& r_force_vector = *mpForceVector;
        // // const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();
        
        // boost::timer force_vector_build_time;
        // if (r_force_vector.size() != system_size)
        //     r_force_vector.resize(system_size, false);
        // r_force_vector = ZeroVector( system_size );
        // p_builder_and_solver->BuildRHS(p_scheme,r_model_part,r_force_vector);
        
        // if (BaseType::GetEchoLevel() > 0 && rank == 0)
        // {
        //     std::cout << "force_vector_build_time : " << force_vector_build_time.elapsed() << std::endl;
        // }

        // // initialize the modal matrix
        // auto& r_modal_matrix = *mpModalMatrix;
        // const std::size_t n_modes = r_process_info[EIGENVALUE_VECTOR].size();
        // if( r_modal_matrix.size1() != system_size || r_modal_matrix.size2() != n_modes )
        //     r_modal_matrix.resize( system_size, n_modes, false );
        // r_modal_matrix = ZeroMatrix( system_size, n_modes );

        // boost::timer modal_matrix_build_time;
        // for( std::size_t i = 0; i < n_modes; ++i )
        // {
        //     for( auto& node : r_model_part.Nodes() )
        //     {
        //         ModelPart::NodeType::DofsContainerType node_dofs = node.GetDofs();
        //         const std::size_t n_node_dofs = node_dofs.size();
        //         const Matrix& r_node_eigenvectors = node.GetValue(EIGENVECTOR_MATRIX);

        //         if( node_dofs.IsSorted() == false )
        //         {
        //             node_dofs.Sort();
        //         }

        //         for( std::size_t j = 0; j < n_node_dofs; ++j )
        //         {
        //             const auto it_dof = std::begin(node_dofs) + j;
        //             r_modal_matrix(it_dof->EquationId(), i) = r_node_eigenvectors(i, j);
        //         }
        //     }
        // }

        // if (BaseType::GetEchoLevel() > 0 && rank == 0)
        // {
        //     std::cout << "modal_matrix_build_time : " << modal_matrix_build_time.elapsed() << std::endl;
        // }

        // // get the damping coefficients if they exist
        // for( auto& property : r_model_part.PropertiesArray() )
        // {
        //     if( property->Has(SYSTEM_DAMPING_RATIO) )
        //     {
        //         mSystemDamping = property->GetValue(SYSTEM_DAMPING_RATIO);
        //     }
            
        //     if( property->Has(RAYLEIGH_ALPHA) && property->Has(RAYLEIGH_BETA) )
        //     {
        //         mRayleighAlpha = property->GetValue(RAYLEIGH_ALPHA);
        //         mRayleighBeta = property->GetValue(RAYLEIGH_BETA);
        //     }
        // }

        // // compute the effective material damping if required
        // if( mUseMaterialDamping )
        // {
        //     // throw an error, if no submodelparts are present
        //     KRATOS_ERROR_IF(r_model_part.NumberOfSubModelParts() < 1) << "No submodelparts detected!" << std::endl;
            
        //     //initialize all required variables
        //     r_model_part.GetProcessInfo()[BUILD_LEVEL] = 2;
        //     mMaterialDampingRatios = ZeroVector( n_modes );
            
        //     //initialize dummy vectors
        //     auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        //     auto pb = SparseSpaceType::CreateEmptyVectorPointer();
        //     auto& rDx = *pDx;
        //     auto& rb = *pb;
        //     SparseSpaceType::Resize(rDx,system_size);
        //     SparseSpaceType::Set(rDx,0.0);
        //     SparseSpaceType::Resize(rb,system_size);
        //     SparseSpaceType::Set(rb,0.0);

        //     //loop over all modes and initialize the material damping ratio per mode
        //     boost::timer material_damping_build_time;
        
        //     for( std::size_t i = 0; i < n_modes; ++i )
        //     {
        //         double up = 0.0;
        //         double down = 0.0;
        //         auto modal_vector = column( r_modal_matrix, i );
        //         for( auto& sub_model_part : r_model_part.SubModelParts() )
        //         {
        //             double damping_coefficient = 0.0;
        //             for( auto& property : sub_model_part.PropertiesArray() )
        //             {
        //                 if( property->Has(SYSTEM_DAMPING_RATIO) )
        //                 {
        //                     damping_coefficient = property->GetValue(SYSTEM_DAMPING_RATIO);
        //                 }
        //             }
                    
        //             //initialize the submodelpart stiffness matrix
        //             auto temp_stiffness_matrix = SparseSpaceType::CreateEmptyMatrixPointer();
        //             p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, 
        //                 temp_stiffness_matrix,
        //                 pDx,
        //                 pb,
        //                 r_model_part.Elements(),
        //                 r_model_part.Conditions(),
        //                 r_model_part.GetProcessInfo());

        //             //build stiffness matrix for submodelpart material
        //             p_builder_and_solver->BuildLHS(p_scheme, sub_model_part, *temp_stiffness_matrix);

        //             //compute strain energy of the submodelpart and the effective damping ratio
        //             double strain_energy = 0.5 * inner_prod( prod(modal_vector, *temp_stiffness_matrix), modal_vector );
        //             down += strain_energy;
        //             up += damping_coefficient * strain_energy;
        //         }
        //         KRATOS_ERROR_IF(down < std::numeric_limits<double>::epsilon()) << "No valid effective "
        //             << "material damping ratio could be computed. Are all elements to be damped available "
        //             << "in the submodelparts? Are the modal vectors available?" << std::endl;
                
        //         mMaterialDampingRatios(i) = up / down;
        //     }

        //     if (BaseType::GetEchoLevel() > 0 && rank == 0)
        //     {
        //         std::cout << "modal_matrix_build_time : " << material_damping_build_time.elapsed() << std::endl;
        //         KRATOS_WATCH(mMaterialDampingRatios)
        //     }
        // }
        
        KRATOS_CATCH("")
    }

    double Solve() override
    {
        KRATOS_TRY

        // auto& r_model_part = BaseType::GetModelPart();

        // operations to be done once
        if (this->GetIsInitialized() == false)
        {
            Initialize();
            this->SetIsInitialized(true);
        }

        this->InitializeSolutionStep();

        // auto& r_process_info = r_model_part.GetProcessInfo();
        // double excitation_frequency = r_process_info[TIME];

        // // get eigenvalues
        // DenseVectorType eigenvalues = r_process_info[EIGENVALUE_VECTOR];
        // const std::size_t n_modes = eigenvalues.size();

        // // DenseMatrixType eigenvectors;
        // const std::size_t n_dofs = this->pGetBuilderAndSolver()->GetEquationSystemSize();
        
        // auto& f = *mpForceVector;

        // ComplexType mode_weight;
        // ComplexVectorType modal_displacement;
        // modal_displacement.resize(n_dofs, false);
        // modal_displacement = ZeroVector( n_dofs );

        // double modal_damping = 0.0;

        // for( std::size_t i = 0; i < n_modes; ++i )
        // {
        //     KRATOS_ERROR_IF(eigenvalues[i] < std::numeric_limits<double>::epsilon()) << "No valid eigenvalue "
        //             << "for mode " << i << std::endl;
        //     modal_damping = mSystemDamping + mRayleighAlpha / (2 * eigenvalues[i]) + mRayleighBeta * eigenvalues[i] / 2;
            
        //     if( mUseMaterialDamping )
        //     {
        //         modal_damping += mMaterialDampingRatios[i];
        //     }

        //     auto& r_modal_matrix = *mpModalMatrix;

        //     DenseVectorType modal_vector(n_dofs);
        //     TDenseSpace::GetColumn(i, r_modal_matrix, modal_vector);

        //     ComplexType factor( eigenvalues[i] - std::pow( excitation_frequency, 2.0 ), 2 * modal_damping * std::sqrt(eigenvalues[i]) * excitation_frequency );
        //     KRATOS_ERROR_IF( std::abs(factor) < std::numeric_limits<double>::epsilon() ) << "No valid modal weight" << std::endl;
        //     mode_weight = inner_prod( modal_vector, f ) / factor;

        //     // compute the modal displacement as a superposition of modal_weight * eigenvector
        //     for( auto& node : r_model_part.Nodes() )
        //     {
        //         auto& node_dofs = node.GetDofs();
        //         const std::size_t n_node_dofs = node_dofs.size();
        //         const Matrix& r_node_eigenvectors = node.GetValue(EIGENVECTOR_MATRIX);

        //         if (node_dofs.IsSorted() == false)
        //         {
        //             node_dofs.Sort();
        //         }

        //         for (std::size_t j = 0; j < n_node_dofs; j++)
        //         {
        //             auto it_dof = std::begin(node_dofs) + j;
        //             modal_displacement[it_dof->EquationId()] = modal_displacement[it_dof->EquationId()] + mode_weight * r_node_eigenvectors(i,j);
        //         }
        //     }
        // }

        // this->AssignVariables(modal_displacement);
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
    void AssignVariables(ComplexVectorType& rModalDisplacement, int step=0)
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
                    it_dof->GetSolutionStepValue(step) = std::abs(rModalDisplacement(it_dof->EquationId()));
                    //phase angle
                    it_dof->GetSolutionStepReactionValue(step) = std::abs(std::arg(rModalDisplacement(it_dof->EquationId())));
                }
                else
                {
                    it_dof->GetSolutionStepValue(step) = 0.0;
                }
            }
        }
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

