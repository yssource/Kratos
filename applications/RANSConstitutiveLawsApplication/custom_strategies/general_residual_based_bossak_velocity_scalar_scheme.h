//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED

// System includes

// Project includes
#include "includes/model_part.h"
#include "solving_strategies/schemes/residual_based_bossak_velocity_scheme.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_constitutive_laws_application_variables.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakVelocityScalarScheme
    : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakVelocityScalarScheme);

    typedef Node<3> NodeType;

    typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::SystemMatrixType SystemMatrixType;

    typedef typename BaseType::SystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    /// Constructor.

    ResidualBasedBossakVelocityScalarScheme(const double AlphaBossak,
                                            Variable<double> const * pScalarVariable,
                                            Variable<double> const * pScalarRateVariable,
                                            Variable<double> const * pRelaxedScalarRateVariable)
        : ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(AlphaBossak, {}, {pScalarVariable}, {pScalarRateVariable}, {}, {}, {}),
        mpScalarRateVariable(pScalarRateVariable),
        mpRelaxedScalarRateVariable(pRelaxedScalarRateVariable)
    {
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::Update(rModelPart, rDofSet, rA, rDx, rb);

        // Updating the auxiliary variables
        const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int iNode = 0; iNode < number_of_nodes; ++iNode)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + iNode);
            const double scalar_rate_dot_old = r_node.FastGetSolutionStepValue(
                *mpScalarRateVariable, 1);
            const double scalar_rate_dot = r_node.FastGetSolutionStepValue(
                *mpScalarRateVariable, 0);

            r_node.FastGetSolutionStepValue(*mpRelaxedScalarRateVariable) =
                this->mAlphaBossak * scalar_rate_dot_old + (1.0 - this->mAlphaBossak) * scalar_rate_dot;
        }

        KRATOS_CATCH("");
    }

private:
    Variable<double> const * mpScalarRateVariable;
    Variable<double> const * mpRelaxedScalarRateVariable;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED defined