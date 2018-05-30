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

            std::vector<double> dummy_vec(5);
            comm_to_test.SumAll(dummy_vec);

        }

#ifdef KRATOS_DEBUG
        KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorExtensionTest2, KratosMappingAppFastSuite)
        {
            ModelPart dummy_model_part("dummy");

            MPICommunicator comm_to_test(&dummy_model_part.GetNodalSolutionStepVariablesList());

            std::vector<double> dummy_vec(5);
            KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(comm_to_test.SumAll(dummy_vec),
                "Error: The ModelPart named : \"Random\" was not found as SubModelPart of : \"Main\". The total input string was \"Main.Random\"");

        }
#endif

#endif

    } // namespace Testing
}  // namespace Kratos.
