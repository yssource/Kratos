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

#if !defined(KRATOS_RANS_VARIABLE_UTILS)
#define KRATOS_RANS_VARIABLE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

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

/**
 * @class VariableUtils
 * @ingroup KratosCore
 * @brief This class implements a set of auxiliar, already parallelized, methods to
 * perform some common tasks related with the variable values and fixity.
 * @details The methods are exported to python in order to add this improvements to the python interface
 * @author Riccardo Rossi
 * @author Ruben Zorrilla
 * @author Vicente Mataix Ferrandiz
 */
class RansVariableUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// We create the Pointer related to RansVariableUtils
    KRATOS_CLASS_POINTER_DEFINITION(RansVariableUtils);

    template <class TVarType>
    void AddToHistoricalNodeScalarVariable(const TVarType& rDestinationVar,
                                           const TVarType& rVarA,
                                           const TVarType& rVarB,
                                           ModelPart& rModelPart,
                                           const unsigned int rBuffStep = 0) const
    {
        KRATOS_TRY

        const int number_of_nodes =
            rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            auto it_node = rModelPart.GetCommunicator().LocalMesh().NodesBegin() + i;
            const double val_a = it_node->GetSolutionStepValue(rVarA, rBuffStep);
            const double val_b = it_node->GetSolutionStepValue(rVarB, rBuffStep);
            double& val_destination =
                it_node->GetSolutionStepValue(rDestinationVar, rBuffStep);
            val_destination = val_a + val_b;
        }

        KRATOS_CATCH("")
    }

    void FixScalarVariableDofs(const Flags& rFlag,
                               const Variable<double>& rVariable,
                               ModelPart& rModelPart) const
    {
        KRATOS_TRY

        const int number_of_nodes =
            rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            auto& r_node = *(rModelPart.GetCommunicator().LocalMesh().NodesBegin() + i);
            if (r_node.Is(rFlag))
                r_node.Fix(rVariable);
        }

        KRATOS_CATCH("")
    }

    unsigned int GetNumberOfNegativeScalarValueNodes(const ModelPart& rModelPart,
                                                     const Variable<double>& rVariable) const
    {
        const int number_of_nodes = rModelPart.NumberOfNodes();

        unsigned int number_of_negative_nodes = 0;

#pragma omp parallel for reduction(+ : number_of_negative_nodes)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const double value =
                (rModelPart.NodesBegin() + i)->FastGetSolutionStepValue(rVariable);
            if (value < 0.0)
            {
                number_of_negative_nodes++;
            }
        }

        return number_of_negative_nodes;
    }

    double GetMinimumScalarValue(const ModelPart& rModelPart,
                                 const Variable<double>& rVariable) const
    {
        const int number_of_nodes = rModelPart.NumberOfNodes();

        if (number_of_nodes == 0)
            return 0.0;

        double min_value = rModelPart.NodesBegin()->FastGetSolutionStepValue(rVariable);

#pragma omp parallel for reduction(min : min_value)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const double value =
                (rModelPart.NodesBegin() + i)->FastGetSolutionStepValue(rVariable);
            min_value = std::min(min_value, value);
        }

        return min_value;
    }

    double GetMaximumScalarValue(const ModelPart& rModelPart,
                                 const Variable<double>& rVariable) const
    {
        const int number_of_nodes = rModelPart.NumberOfNodes();

        if (number_of_nodes == 0)
            return 0.0;

        double max_value = rModelPart.NodesBegin()->FastGetSolutionStepValue(rVariable);

#pragma omp parallel for reduction(max : max_value)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const double value =
                (rModelPart.NodesBegin() + i)->FastGetSolutionStepValue(rVariable);
            max_value = std::max(max_value, value);
        }

        return max_value;
    }

    double GetScalarVariableIncreaseNormSquare(const ModelPart& rModelPart,
                                               const Variable<double>& rOldVariable,
                                               const Variable<double>& rNewVariable) const
    {
        const int number_of_nodes = rModelPart.NumberOfNodes();

        double increase_norm_square = 0.0;

#pragma omp parallel for reduction(+ : increase_norm_square)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const auto& r_node = *(rModelPart.NodesBegin() + i);
            const double old_value = r_node.FastGetSolutionStepValue(rOldVariable);
            const double new_value = r_node.FastGetSolutionStepValue(rNewVariable);
            increase_norm_square += std::pow(old_value - new_value, 2);
        }

        return increase_norm_square;
    }

    double GetScalarVariableSolutionNormSquare(const ModelPart& rModelPart,
                                               const Variable<double>& rVariable) const
    {
        const int number_of_nodes = rModelPart.NumberOfNodes();

        double solution_norm_square = 0.0;

#pragma omp parallel for reduction(+ : solution_norm_square)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const double solution_value =
                (rModelPart.NodesBegin() + i)->FastGetSolutionStepValue(rVariable);
            solution_norm_square += std::pow(solution_value, 2);
        }

        return solution_norm_square;
    }
}; /* Class RansVariableUtils */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RANS_VARIABLE_UTILS  defined */