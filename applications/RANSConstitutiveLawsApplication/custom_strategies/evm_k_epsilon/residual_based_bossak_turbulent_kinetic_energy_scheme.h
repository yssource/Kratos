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

    typedef typename BaseType::SystemMatrixType SystemMatrixType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::SystemVectorType SystemVectorType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

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

    // void InitializeSolutionStep(
    //     ModelPart& rModelPart,
    //     SystemMatrixType& A,
    //     SystemVectorType& Dx,
    //     SystemVectorType& b
    // ) override
    // {
    //     ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    //     r_current_process_info[RANS_STABILIZATION_MULTIPLIER] = 1.0;
    // }

    // void FinalizeNonLinIteration(
    //     ModelPart& rModelPart,
    //     SystemMatrixType& A,
    //     SystemVectorType& Dx,
    //     SystemVectorType& b
    //     ) override
    // {
    //     ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    //     double& stab_multiplier = r_current_process_info[RANS_STABILIZATION_MULTIPLIER];
    //     const double stab_multiplier_max = r_current_process_info[RANS_STABILIZATION_MULTIPLIER_MAX];

    //     if (mPreviousStabilizationMultiplier == 0.0)
    //     {
    //         stab_multiplier += 1.0;
    //         mPreviousStabilizationMultiplier = stab_multiplier;
    //         mPreviousAggregatedError = CalculationUtilities::WarnIfNegative<NodeType>(rModelPart, TURBULENT_KINETIC_ENERGY, "KScheme");
    //     } else
    //     {
    //         const double current_aggregated_error = CalculationUtilities::WarnIfNegative<NodeType>(rModelPart, TURBULENT_KINETIC_ENERGY, "KScheme");
    //         if (current_aggregated_error != 0.0)
    //         {
    //             const double m = (current_aggregated_error - mPreviousAggregatedError) / (stab_multiplier - mPreviousStabilizationMultiplier);
    //             mPreviousStabilizationMultiplier = stab_multiplier;
    //             mPreviousAggregatedError = current_aggregated_error;
    //             const double c = current_aggregated_error - m * stab_multiplier;
    //             stab_multiplier = std::min(stab_multiplier_max, -c / m);
    //         }
    //     }

    // }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        BaseType::CalculateSystemContributions(pCurrentElement, rLHS_Contribution,
                                               rRHS_Contribution, rEquationId,
                                               rCurrentProcessInfo);

        // const int k = OpenMPUtils::ThisThread();
        // this->mDampingMatrix[k].clear();

        // for (unsigned int i = 0; i < rLHS_Contribution.size1(); ++i)
        // {
        //     for (unsigned int j = 0; j < rLHS_Contribution.size2(); ++j)
        //     {
        //         const double max_value =
        //             std::max(rLHS_Contribution(i, j), rLHS_Contribution(j, i));
        //         this->mDampingMatrix[k](i, j) = -std::max(0.0, max_value);
        //     }
        // }

        // for (unsigned int i = 0; i < rLHS_Contribution.size1(); ++i)
        // {
        //     this->mDampingMatrix[k](i, i) = 0.0;
        //     for (unsigned int j = 0; j < rLHS_Contribution.size2(); ++j)
        //     {
        //         this->mDampingMatrix[k](i, i) += this->mDampingMatrix[k](i, j);
        //     }
        // }

        // double residual = 0.0;
        // pCurrentElement->Calculate(RESIDUAL, residual, rCurrentProcessInfo);

        // noalias(this->mDampingMatrix[k]) = this->mDampingMatrix[k] * residual;

        // Vector U = ZeroVector(rLHS_Contribution.size1());
        // for (unsigned int iNode = 0; iNode < rLHS_Contribution.size1(); ++iNode)
        //     U[iNode] = pCurrentElement->GetGeometry()[iNode].FastGetSolutionStepValue(
        //         TURBULENT_KINETIC_ENERGY);

        // noalias(rLHS_Contribution) += this->mDampingMatrix[k];

        // noalias(rRHS_Contribution) -= prod(this->mDampingMatrix[k], U);

        KRATOS_CATCH("");
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
            const double tke_dot_old =
                r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 1);
            const double tke_dot =
                r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 0);

            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) =
                this->mAlphaBossak * tke_dot_old + (1.0 - this->mAlphaBossak) * tke_dot;
        }

        KRATOS_CATCH("");
    }

    ///@}

private:
    double mPreviousStabilizationMultiplier = 0.0;
    double mPreviousAggregatedError = 0.0;
};

} // namespace Kratos

#endif // KRATOS_RESIDUAL_BASED_BOSSAK_TURBULENT_KINETIC_ENERGY_SCHEME_H_INCLUDED defined