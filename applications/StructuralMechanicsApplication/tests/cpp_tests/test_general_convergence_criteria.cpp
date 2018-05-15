// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Natalia Saiapova
//                   Philipp Bucher
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

#include "custom_strategies/custom_convergencecriterias/general_residual_criteria.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
// #include "includes/mpi_communicator.h"
#endif

/*
Things we will test:
- For both residual- & solutionupdate-based => here we can think if it is really necessary for both, maybe for one of them just one test if they have a common baseclass
    - For relative and absolute convergence
        - For different types of variables:
            - No Input (e.g. what the Displacement and the Residual Criterion do atm)
            - One Array3D + Double Variable (what VelPrCriteria does)
            - Two Double Variables
            - Two Array3D Variables (What the criteria in StrucutralMechanics do atm)
            - ... ?
=> makes at least 2*2*(4+x) tests
*/

namespace Kratos
{
    namespace Testing
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;


        typedef typename SparseSpaceType::DataType TDataType;

        typedef typename SparseSpaceType::MatrixType TSystemMatrixType;
        typedef typename SparseSpaceType::VectorType TSystemVectorType;
        typedef typename SparseSpaceType::MatrixPointerType TSystemMatrixPointerType;
        typedef typename SparseSpaceType::VectorPointerType TSystemVectorPointerType;

        typedef typename LocalSpaceType::MatrixType TLocalMatrixType;
        typedef typename LocalSpaceType::VectorType TLocalVectorType;
        typedef typename LocalSpaceType::MatrixPointerType TLocalMatrixPointerType;
        typedef typename LocalSpaceType::VectorPointerType TLocalVectorPointerType;


        typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
        typedef Kratos::unique_ptr<ConvergenceCriteriaType> ConvergenceCriteriaPointerType;

        typedef GeneralResidualCriteria<SparseSpaceType, LocalSpaceType> GenConvergenceCriteriaType;
        typedef Kratos::unique_ptr<GenConvergenceCriteriaType> GenConvergenceCriteriaPointerType;

        typedef ModelPart::DofsArrayType DofsArrayType;

        typedef Node<3> NodeType;
        typedef Kratos::unique_ptr<NodeType> NodeUniquePointerType;

        // This function initializes the test depending on the passed configuration
        // TODO for now the dofs are hard-coded, this should be made variable
        void SetUpTest(DofsArrayType& rDofSet,
                       TSystemVectorType& rSystemVec_Dx,
                       TSystemVectorType& rSystemVec_b,
                       const std::size_t SystemSize)
        {
            std::vector<NodeUniquePointerType> node_pointer_vector(SystemSize);

            const std::size_t num_dofs_per_node = 7; // TODO hard coded for now, to be made variable

            std::vector< Dof<double>::Pointer > dof_pointer_vector;

            dof_pointer_vector.reserve((SystemSize*num_dofs_per_node));

            // Creating Nodes and adding Dofs
            for (std::size_t i=0; i<SystemSize; ++i)
            {
                node_pointer_vector[i] = Kratos::make_unique<NodeType>();

                node_pointer_vector[i]->AddDof(ROTATION_X, REACTION_MOMENT_X);
                node_pointer_vector[i]->AddDof(ROTATION_Y, REACTION_MOMENT_Y);
                node_pointer_vector[i]->AddDof(ROTATION_Z, REACTION_MOMENT_Z);

                node_pointer_vector[i]->AddDof(VELOCITY_X, REACTION_X);
                node_pointer_vector[i]->AddDof(VELOCITY_Y, REACTION_Y);
                node_pointer_vector[i]->AddDof(VELOCITY_Z, REACTION_Z);

                node_pointer_vector[i]->AddDof(PRESSURE, REACTION_WATER_PRESSURE);

                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(ROTATION_X));
                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(ROTATION_Y));
                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(ROTATION_Z));
                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(VELOCITY_X));
                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(VELOCITY_Y));
                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(VELOCITY_Z));
                dof_pointer_vector.push_back(node_pointer_vector[i]->pGetDof(PRESSURE));
            }

            rDofSet.clear();
            rDofSet.reserve(dof_pointer_vector.size());

            for (std::size_t i=0; i<dof_pointer_vector.size(); ++i)
                dof_pointer_vector[i]->SetEquationId(i);

            for (auto it=dof_pointer_vector.begin(); it!=dof_pointer_vector.end(); it++)
                rDofSet.push_back( it->get() );

            rDofSet.Sort();

            // Initializing the Solution Vector
            SparseSpaceType::Resize(rSystemVec_Dx, SystemSize*num_dofs_per_node);
            SparseSpaceType::SetToZero(rSystemVec_Dx);
            SparseSpaceType::Resize(rSystemVec_b, SystemSize*num_dofs_per_node);
            SparseSpaceType::SetToZero(rSystemVec_b);
        }

        // This function sets the solution on the dofs and the system vectors
        // Usually this would be done by setting up and solving the system
        // Natasha please implement this, think abt sth smart to set an apropriate solution
        void SetSolution(DofsArrayType& rDofSet,
                         TSystemVectorType& rSystemVec_Dx,
                         TSystemVectorType& rSystemVec_b)
        {
            // for (auto& r_dof : rDofSet)
            //     r_dof.GetSolutionStepValue() = 5.123; // TODO so far this is a random value
            // This does not work yet, I have to first figure out why ...

            const std::size_t system_size = SparseSpaceType::Size(rSystemVec_Dx);
            KRATOS_ERROR_IF(system_size != SparseSpaceType::Size(rSystemVec_b))
                << "Sytsem Vector sizes are inconsistent!" << std::endl;

            for (std::size_t i=0; i< system_size; ++i)
            {
                rSystemVec_Dx[i] = -5.147;
                rSystemVec_b[i] = 565.147;
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(GeneralResudialConvergenceCriteriaTest1, KratosStructuralMechanicsFastSuite)
        {
            DofsArrayType dofs_array;
            const std::size_t system_size = 10;

            // Natasha even though I am initializing both b and Dx the implementation of the criteria should
            // for now only consider b!
            // In a second step we will make it more general that it can be used with either b or Dx
            TSystemMatrixPointerType p_dummy_system_matrix = SparseSpaceType::CreateEmptyMatrixPointer();
            TSystemVectorPointerType p_system_vector_b = SparseSpaceType::CreateEmptyVectorPointer();
            TSystemVectorPointerType p_system_vector_Dx = SparseSpaceType::CreateEmptyVectorPointer();

            TSystemMatrixType& r_system_matrix = *p_dummy_system_matrix;
            TSystemVectorType& r_system_vector_b = *p_system_vector_b;
            TSystemVectorType& r_system_vector_Dx = *p_system_vector_Dx;

            SetUpTest(dofs_array, r_system_vector_Dx, r_system_vector_b, system_size);

            const TDataType NewRatioTolerance = 1e-7;
            const TDataType AlwaysConvergedNorm = 1e-5;

            ConvergenceCriteriaPointerType p_conv_crit = Kratos::make_unique<GenConvergenceCriteriaType>(
                NewRatioTolerance, AlwaysConvergedNorm);

            ModelPart dummy_model_part("dummy");

            const std::size_t num_solution_steps = 5; // Natasha this should be used by "SetSolution"
            // in order to set a solution that will make the problem converge

            p_conv_crit->Initialize(dummy_model_part);

            bool is_converged = false;

            std::size_t counter = 0;

            while(!is_converged)
            {
                SetSolution(dofs_array, r_system_vector_Dx, r_system_vector_b);
                p_conv_crit->InitializeSolutionStep(dummy_model_part,
                                                    dofs_array,
                                                    r_system_matrix,
                                                    r_system_vector_Dx,
                                                    r_system_vector_b);

                is_converged = p_conv_crit->PostCriteria(dummy_model_part,
                                                    dofs_array,
                                                    r_system_matrix,
                                                    r_system_vector_Dx,
                                                    r_system_vector_b);

                counter += 1;

                if (counter > 10) // TODO this is an intermediate solution until the criteria is implemented
                    is_converged = true;
            }
        }

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorExtensionTest, KratosStructuralMechanicsFastSuite)
        {
            ModelPart dummy_model_part("dummy");

            const auto& dummy_var_list = dummy_model_part.GetNodalSolutionStepVariablesList();

            // MPICommunicator comm_to_test(*());

            // comm_to_test.SumAll(...);

        }
#endif

    } // namespace Testing
}  // namespace Kratos.
