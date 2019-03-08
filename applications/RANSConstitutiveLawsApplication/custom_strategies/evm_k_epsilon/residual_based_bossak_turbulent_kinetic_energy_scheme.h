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

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_KINETIC_ENERGY_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_KINETIC_ENERGY_SCHEME_H_INCLUDED

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
class ResidualBasedBossakTurbulentKineticEnergyScheme
    : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>
{
    typedef Node<3> NodeType;

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
            NodeType& r_node = mpElement->GetGeometry()[NodeId];
            rVector[0] = MakeIndirectScalar(r_node, TURBULENT_KINETIC_ENERGY, Step);
        }

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step,
                                        ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            NodeType& r_node = mpElement->GetGeometry()[NodeId];
            rVector[0] = MakeIndirectScalar(r_node, TURBULENT_KINETIC_ENERGY_RATE, Step);
        }

        void GetFirstDerivativesDofsVector(std::size_t NodeId,
                                           std::vector<Dof<double>::Pointer>& rVector,
                                           ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            NodeType& r_node = mpElement->GetGeometry()[NodeId];
            rVector[0] = r_node.pGetDof(TURBULENT_KINETIC_ENERGY);
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                          ProcessInfo& rCurrentProcessInfo) const override
        {
            rVariables.resize(1);
            rVariables[0] = &TURBULENT_KINETIC_ENERGY;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                           ProcessInfo& rCurrentProcessInfo) const override
        {
            rVariables.resize(1);
            rVariables[0] = &TURBULENT_KINETIC_ENERGY_RATE;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakTurbulentKineticEnergyScheme);

    typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Constructor.

    ResidualBasedBossakTurbulentKineticEnergyScheme(const double AlphaBossak)
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
            Element& r_element = *(rModelPart.ElementsBegin() + i);
            r_element.SetValue(DERIVATIVES_EXTENSION,
                               Kratos::make_shared<ElementDerivativesExtension>(&r_element));
        }

        KRATOS_INFO("KScheme") << "Initialized.\n";

        KRATOS_CATCH("");
    }
    ///@}
};

} // namespace Kratos

#endif // KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_KINETIC_ENERGY_SCHEME_H_INCLUDED defined