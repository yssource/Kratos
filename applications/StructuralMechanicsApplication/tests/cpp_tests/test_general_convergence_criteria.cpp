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


namespace Kratos
{
    namespace Testing
    {

        KRATOS_TEST_CASE_IN_SUITE(GeneralResudialConvergenceCriteriaTest1, KratosStructuralMechanicsFastSuite)
        {
            // ModelPart this_model_part("Main");
            // this_model_part.SetBufferSize(2);

            // PrismNeighboursProcessCreateModelPart(this_model_part, 3);

            // PrismNeighboursProcess prism_neighbours_process = PrismNeighboursProcess(this_model_part);
            // prism_neighbours_process.Execute();

            // auto pneigh = (this_model_part.Elements().begin())->GetValue(NEIGHBOUR_NODES);
            // KRATOS_CHECK_EQUAL(pneigh[0].Id(), 7);
            // KRATOS_CHECK_EQUAL(pneigh[1].Id(), 9);
            // KRATOS_CHECK_EQUAL(pneigh[2].Id(), 8);
        }

    } // namespace Testing
}  // namespace Kratos.
