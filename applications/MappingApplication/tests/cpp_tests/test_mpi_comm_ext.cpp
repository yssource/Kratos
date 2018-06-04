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


/*
References:
- applications/mpi_search_application/tests/test_global_pointer_mpi.cpp
- kratos/includes/checks.h

*/

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "includes/mpi_communicator.h"
#endif

namespace Kratos
{
    namespace Testing
    {

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorExtensionTest, KratosMappingAppFastSuite)
        {
            ModelPart dummy_model_part("dummy");

            MPICommunicator comm_to_test(&dummy_model_part.GetNodalSolutionStepVariablesList());
            int my_pid = comm_to_test.MyPID();
            int n_ranks = comm_to_test.TotalProcesses();
            double result_sum = n_ranks * (n_ranks - 1) / 2;

            std::vector<double> dummy_vec(5, (double) my_pid);
            comm_to_test.SumAll(dummy_vec);

            std::for_each(dummy_vec.begin(), dummy_vec.end(), [result_sum](double el) {
                KRATOS_CHECK_DOUBLE_EQUAL(el, result_sum);
            });

        }

#ifdef KRATOS_DEBUG
        KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorExtensionTest2, KratosMappingAppFastSuite)
        {
            ModelPart dummy_model_part("dummy");

            MPICommunicator comm_to_test(&dummy_model_part.GetNodalSolutionStepVariablesList());
            int my_pid = comm_to_test.MyPID();
            std::vector<double> dummy_vec;

            if (my_pid % 2 == 0) {
                dummy_vec.resize(6);
            } else {
                dummy_vec.resize(5);
            }

            KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(comm_to_test.SumAll(dummy_vec),
                "Error: The ModelPart named : \"Random\" was not found as SubModelPart of : \"Main\". The total input string was \"Main.Random\"");

        }
#endif

#endif

    } // namespace Testing
}  // namespace Kratos.
