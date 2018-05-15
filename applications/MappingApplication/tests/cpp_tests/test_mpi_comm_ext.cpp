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

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "includes/mpi_communicator.h"
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

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorExtensionTest, KratosMappingAppFastSuite)
        {
            ModelPart dummy_model_part("dummy");

            MPICommunicator comm_to_test(&dummy_model_part.GetNodalSolutionStepVariablesList());

            std::vector<double> dummy_vec(5);
            comm_to_test.SumAll(dummy_vec);

        }
#endif

    } // namespace Testing
}  // namespace Kratos.
