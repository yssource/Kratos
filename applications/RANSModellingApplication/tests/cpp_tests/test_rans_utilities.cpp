//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <random>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/ublas_interface.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RansAssembleElementMatrix, RansCalculationUtilities)
{
    const int number_of_nodes = 3;

    const int elem_1_dimension = 4;
    const int elem_1_size = elem_1_dimension * number_of_nodes;

    const int elem_2_dimension = 2;
    const int elem_2_size = elem_2_dimension * number_of_nodes;

    const int elem_3_dimension = 3;
    const int elem_3_size = elem_3_dimension * number_of_nodes;

    const int total_block_size = (elem_1_dimension + elem_2_dimension + elem_3_dimension);
    const int total_size = number_of_nodes * total_block_size;

    Matrix assembled_matrix(total_size, total_size);

    Matrix elem_1_1(elem_1_size, elem_1_size);
    Matrix elem_1_2(elem_2_size, elem_1_size);
    Matrix elem_1_3(elem_3_size, elem_1_size);

    Matrix elem_2_1(elem_1_size, elem_2_size);
    Matrix elem_2_2(elem_2_size, elem_2_size);
    Matrix elem_2_3(elem_3_size, elem_2_size);

    Matrix elem_3_1(elem_1_size, elem_3_size);
    Matrix elem_3_2(elem_2_size, elem_3_size);
    Matrix elem_3_3(elem_3_size, elem_3_size);

    int local_value = 1;

    for (int row_node = 0; row_node < number_of_nodes; ++row_node)
    {
        for (int column_node = 0; column_node < number_of_nodes; ++column_node)
        {
            for (int elem_1_row_dim = 0; elem_1_row_dim < elem_1_dimension; ++elem_1_row_dim)
            {
                for (int elem_1_column_dim = 0;
                     elem_1_column_dim < elem_1_dimension; ++elem_1_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_row_dim,
                                     column_node * total_block_size + elem_1_column_dim) =
                        local_value;
                    elem_1_1(row_node * elem_1_dimension + elem_1_row_dim,
                             column_node * elem_1_dimension + elem_1_column_dim) = local_value;
                }
            }
            for (int elem_2_row_dim = 0; elem_2_row_dim < elem_2_dimension; ++elem_2_row_dim)
            {
                for (int elem_1_column_dim = 0;
                     elem_1_column_dim < elem_1_dimension; ++elem_1_column_dim)
                {
                    local_value++;
                    assembled_matrix(
                        row_node * total_block_size + elem_1_dimension + elem_2_row_dim,
                        column_node * total_block_size + elem_1_column_dim) = local_value;
                    elem_1_1(row_node * elem_2_dimension + elem_2_row_dim,
                             column_node * elem_1_dimension + elem_1_column_dim) = local_value;
                }
            }
            for (int elem_3_row_dim = 0; elem_3_row_dim < elem_3_dimension; ++elem_3_row_dim)
            {
                for (int elem_1_column_dim = 0;
                     elem_1_column_dim < elem_1_dimension; ++elem_1_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_dimension +
                                         elem_2_dimension + elem_3_row_dim,
                                     column_node * total_block_size + elem_1_column_dim) =
                        local_value;
                    elem_1_1(row_node * elem_3_dimension + elem_3_row_dim,
                             column_node * elem_1_dimension + elem_1_column_dim) = local_value;
                }
            }
            for (int elem_1_row_dim = 0; elem_1_row_dim < elem_1_dimension; ++elem_1_row_dim)
            {
                for (int elem_2_column_dim = 0;
                     elem_2_column_dim < elem_2_dimension; ++elem_2_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_row_dim,
                                     column_node * total_block_size + elem_1_dimension +
                                         elem_2_column_dim) = local_value;
                    elem_1_1(row_node * elem_1_dimension + elem_1_row_dim,
                             column_node * elem_2_dimension + elem_2_column_dim) = local_value;
                }
            }
            for (int elem_2_row_dim = 0; elem_2_row_dim < elem_2_dimension; ++elem_2_row_dim)
            {
                for (int elem_2_column_dim = 0;
                     elem_2_column_dim < elem_2_dimension; ++elem_2_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_dimension + elem_2_row_dim,
                                     column_node * total_block_size + elem_1_dimension +
                                         elem_2_column_dim) = local_value;
                    elem_1_1(row_node * elem_2_dimension + elem_2_row_dim,
                             column_node * elem_2_dimension + elem_2_column_dim) = local_value;
                }
            }
            for (int elem_3_row_dim = 0; elem_3_row_dim < elem_3_dimension; ++elem_3_row_dim)
            {
                for (int elem_2_column_dim = 0;
                     elem_2_column_dim < elem_2_dimension; ++elem_2_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_dimension +
                                         elem_2_dimension + elem_3_row_dim,
                                     column_node * total_block_size + elem_1_dimension +
                                         elem_2_column_dim) = local_value;
                    elem_1_1(row_node * elem_3_dimension + elem_3_row_dim,
                             column_node * elem_2_dimension + elem_2_column_dim) = local_value;
                }
            }
            for (int elem_1_row_dim = 0; elem_1_row_dim < elem_1_dimension; ++elem_1_row_dim)
            {
                for (int elem_3_column_dim = 0;
                     elem_3_column_dim < elem_3_dimension; ++elem_3_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_row_dim,
                                     column_node * total_block_size + elem_1_dimension +
                                         elem_2_dimension + elem_3_column_dim) = local_value;
                    elem_1_1(row_node * elem_1_dimension + elem_1_row_dim,
                             column_node * elem_3_dimension + elem_3_column_dim) = local_value;
                }
            }
            for (int elem_2_row_dim = 0; elem_2_row_dim < elem_2_dimension; ++elem_2_row_dim)
            {
                for (int elem_3_column_dim = 0;
                     elem_3_column_dim < elem_3_dimension; ++elem_3_column_dim)
                {
                    local_value++;
                    assembled_matrix(
                        row_node * total_block_size + elem_1_dimension + elem_2_row_dim,
                        column_node * total_block_size + elem_1_dimension +
                            elem_2_dimension + elem_3_column_dim) = local_value;
                    elem_1_1(row_node * elem_2_dimension + elem_2_row_dim,
                             column_node * elem_3_dimension + elem_3_column_dim) = local_value;
                }
            }
            for (int elem_3_row_dim = 0; elem_3_row_dim < elem_3_dimension; ++elem_3_row_dim)
            {
                for (int elem_3_column_dim = 0;
                     elem_3_column_dim < elem_3_dimension; ++elem_3_column_dim)
                {
                    local_value++;
                    assembled_matrix(row_node * total_block_size + elem_1_dimension +
                                         elem_2_dimension + elem_3_row_dim,
                                     column_node * total_block_size + elem_1_dimension +
                                         elem_2_dimension + elem_3_column_dim) = local_value;
                    elem_1_1(row_node * elem_3_dimension + elem_3_row_dim,
                             column_node * elem_3_dimension + elem_3_column_dim) = local_value;
                }
            }
        }
    }

    RansCalculationUtilities rans_calculation_utilities;

    Matrix check_assembled_matrix(total_size, total_size);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_1_1, number_of_nodes, 0, 0);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_1_2, number_of_nodes, elem_1_dimension, 0);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_1_3, number_of_nodes,
        elem_1_dimension + elem_2_dimension, 0);

    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_2_1, number_of_nodes, 0, elem_1_dimension);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_2_2, number_of_nodes, elem_1_dimension, elem_1_dimension);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_2_3, number_of_nodes,
        elem_1_dimension + elem_2_dimension, elem_1_dimension);

    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_3_1, number_of_nodes, 0, elem_1_dimension + elem_2_dimension);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_3_2, number_of_nodes, elem_1_dimension,
        elem_1_dimension + elem_2_dimension);
    rans_calculation_utilities.AssembleElementMatrix(
        check_assembled_matrix, elem_3_3, number_of_nodes,
        elem_1_dimension + elem_2_dimension, elem_1_dimension + elem_2_dimension);

    RansModellingApplicationTestUtilities::IsMatricesSame(
        check_assembled_matrix, assembled_matrix, std::numeric_limits<double>::epsilon());
}
} // namespace Testing
} // namespace Kratos