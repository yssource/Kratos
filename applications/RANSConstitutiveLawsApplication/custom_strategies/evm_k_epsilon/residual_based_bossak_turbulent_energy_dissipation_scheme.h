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
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_ENERGY_DISSIPATION_RATE_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_ENERGY_DISSIPATION_RATE_SCHEME_H_INCLUDED

// System includes

// Project includes
#include "includes/model_part.h"
#include "solving_strategies/schemes/residual_based_bossak_velocity_scheme.h"
#include "utilities/derivatives_extension.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/calculation_utilities.h"
#include "rans_constitutive_laws_application_variables.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakTurbulentEnergyDissipationRateScheme
    : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>
{
    class ElementDerivativesExtension : public DerivativesExtension
    {
        Element* mpElement;

    public:
        explicit ElementDerivativesExtension(Element* pElement)
            : mpElement(pElement)
        {
        }

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step,
                                       ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            Node<3>& r_node = mpElement->GetGeometry()[NodeId];
            rVector[0] =
                MakeIndirectScalar(r_node, TURBULENT_ENERGY_DISSIPATION_RATE, Step);
        }

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step,
                                        ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            Node<3>& r_node = mpElement->GetGeometry()[NodeId];
            rVector[0] =
                MakeIndirectScalar(r_node, TURBULENT_ENERGY_DISSIPATION_RATE_2, Step);
        }

        void GetFirstDerivativesDofsVector(std::size_t NodeId,
                                           std::vector<Dof<double>::Pointer>& rVector,
                                           ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            Node<3>& r_node = mpElement->GetGeometry()[NodeId];
            rVector[0] = r_node.pGetDof(TURBULENT_ENERGY_DISSIPATION_RATE);
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                          ProcessInfo& rCurrentProcessInfo) const override
        {
            rVariables.resize(1);
            rVariables[0] = &TURBULENT_ENERGY_DISSIPATION_RATE;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                           ProcessInfo& rCurrentProcessInfo) const override
        {
            rVariables.resize(1);
            rVariables[0] = &TURBULENT_ENERGY_DISSIPATION_RATE_2;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakTurbulentEnergyDissipationRateScheme);

    typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Constructor.

    ResidualBasedBossakTurbulentEnergyDissipationRateScheme(const double AlphaBossak)
        : ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(AlphaBossak, true, false)
    {
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        const int number_of_elements = rModelPart.NumberOfElements();

#pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++)
        {
            Element& r_current_element = *(rModelPart.ElementsBegin() + i);
            r_current_element.SetValue(
                DERIVATIVES_EXTENSION,
                Kratos::make_shared<ElementDerivativesExtension>(&r_current_element));
        }

        KRATOS_INFO("EpsilonScheme") << "Initialized.\n";

        KRATOS_CATCH("");
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                 typename BaseType::SystemMatrixType& A,
                                 typename BaseType::SystemVectorType& Dx,
                                 typename BaseType::SystemVectorType& b) override
    {
        KRATOS_TRY

        // LowerBound(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 1e-15);

        KRATOS_CATCH("")
    }

    ///@}
};

} // namespace Kratos

#endif // KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_ENERGY_DISSIPATION_RATE_SCHEME_H_INCLUDED defined