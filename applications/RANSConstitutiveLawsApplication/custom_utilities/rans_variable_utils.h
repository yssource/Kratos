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

#if !defined(KRATOS_RANS_VARIABLE_UTILS )
#define  KRATOS_RANS_VARIABLE_UTILS

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

    /// The nodes container
    typedef ModelPart::NodesContainerType NodesContainerType;

    template< class TVarType >
    void AddToHistoricalNodeScalarVariable(
        const TVarType& rDestinationVar,
        const TVarType& rVarA,
        const TVarType& rVarB,
        const ModelPart& rModelPart,
        const unsigned int rBuffStep = 0
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rModelPart.GetCommunicator().LocalMesh().NumberOfNodes()); ++k) {
            const auto it_node = rModelPart.GetCommunicator().LocalMesh().NodesBegin() + k;
            const double val_a = it_node->GetSolutionStepValue(rVarA, rBuffStep);
            const double val_b = it_node->GetSolutionStepValue(rVarB, rBuffStep);
            double& val_destination = it_node->GetSolutionStepValue(rDestinationVar, rBuffStep);
            val_destination = val_a + val_b;
        }

        KRATOS_CATCH("")
    }

}; /* Class RansVariableUtils */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RANS_VARIABLE_UTILS  defined */