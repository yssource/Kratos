//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_STABILIZED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED)
#define KRATOS_STABILIZED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"
#include "custom_utilities/adjoint_artificial_diffusion.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for dynamic adjoint equations, using Bossak time integration.
/**
 */
template <class TSparseSpace, class TDenseSpace>
class StabilizedAdjointBossakScheme : public ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(StabilizedAdjointBossakScheme);

    typedef ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::SystemVectorType SystemVectorType;
    typedef typename BaseType::SystemMatrixType SystemMatrixType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    StabilizedAdjointBossakScheme(Parameters& rParameters, Parameters& rParametersBossak, AdjointResponseFunction::Pointer pResponseFunction):
        ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace>(rParametersBossak, pResponseFunction),
        mAdjointArtificialDiffusion()
    {
        KRATOS_TRY;

        Parameters default_parameters(R"({
            "scheme_type"                      : "stabilized_bossak",
            "overall_stabilization_coefficient": 0.0,
            "stabilization_source_coefficient" : 1.0,
            "calculate_matrix_energies"        : true,
            "stabilization_settings"           : {},
            "bossak_scheme_settings"           : {}
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mAdjointArtificialDiffusion.SetArtificialDiffusionParameters(rParameters);

        mOverallDiffusionCoefficient = rParameters["overall_stabilization_coefficient"].GetDouble();
        mStabilizationSourceCoefficient = rParameters["stabilization_source_coefficient"].GetDouble();
        mIsMatrixEnergiesCalculated = rParameters["calculate_matrix_energies"].GetBool();

        KRATOS_ERROR_IF(mOverallDiffusionCoefficient < 0.0)<<"Invalid overall stabilization coefficient. \"overall_stabilization_coefficient\" >= 0.0 [ "<<mOverallDiffusionCoefficient<< " < 0.0 ]";
        KRATOS_ERROR_IF(mStabilizationSourceCoefficient < 0.0)<<"Invalid stabilization source coefficient. \"stabilization_source_coefficient\" >= 0.0 [ "<<mStabilizationSourceCoefficient<< " < 0.0 ]";

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

     void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        #pragma omp parallel for
        for (int i = 0;i< static_cast<int>(rModelPart.Elements().size()); ++i)
        {
            auto ie = rModelPart.ElementsBegin() + i;
            ie->SetValue(ARTIFICIAL_DIFFUSION, 0.0);
        }
        // Allocate auxiliary memory.
        int num_threads = OpenMPUtils::GetNumThreads();
        mArtificialDiffusionMatrix.resize(num_threads);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        // Check and resize rLHS and rRHS
        this->CheckAndResizeLocalSystem(pCurrentElement, rLHS_Contribution, rRHS_Contribution);

        // Contribution from variable gradients
        this->CalculateGradientContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contribution from first derivative gradients
        this->CalculateFirstDerivativeContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contribution from second derivative gradients
        this->CalculateSecondDerivativeContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contributions from the previos time step
        this->CalculatePreviousTimeStepContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        LocalSystemMatrixType matrix_29;
        pCurrentElement->Calculate(VMS_ADJOINT_ENERGY_GENERATION_RATE_MATRIX,  matrix_29, rCurrentProcessInfo);

        pCurrentElement->SetValue(STABILIZATION_ANALYSIS_MATRIX_29, matrix_29);

        // Calculate artificial diffusion
        this->CalculateArtificialDiffusionContribution(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Make the local contribution residual
        this->CalculateResidualLocalContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        pCurrentElement->SetValue(STABILIZATION_ANALYSIS_MATRIX_24, rLHS_Contribution);

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(
                        ModelPart& rModelPart,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        if (mIsMatrixEnergiesCalculated)
            {
            ModelPart::ElementsContainerType& r_elements = rModelPart.Elements();
            LocalSystemMatrixType dummy_matrix;

            #pragma omp parallel for
            for (int i=0; i < static_cast<int>(r_elements.size()); ++i)
            {
                auto r_element = r_elements.begin() + i;
                r_element->Calculate(STABILIZATION_ANALYSIS_MATRICES, dummy_matrix, rModelPart.GetProcessInfo());
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateArtificialDiffusionContribution(Element::Pointer pCurrentElement,
                                                  LocalSystemMatrixType& rLHS_Contribution,
                                                  LocalSystemVectorType& rRHS_Contribution,
                                                  ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const Geometry< Node<3> >& r_geometry = pCurrentElement->GetGeometry();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();

        double gauss_integration_weight = 0.0;

        if (domain_size == 2)
            gauss_integration_weight = r_geometry.Area();
        else if (domain_size == 3)
            gauss_integration_weight = r_geometry.Volume();

        double artificial_diffusion;

        int k = OpenMPUtils::ThisThread();
        auto& r_artificial_diffusion_matrix = mArtificialDiffusionMatrix[k];

        artificial_diffusion = mAdjointArtificialDiffusion.CalculateArtificialDiffusion(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
        artificial_diffusion *= mOverallDiffusionCoefficient;
        pCurrentElement->SetValue(ARTIFICIAL_DIFFUSION, artificial_diffusion);
        pCurrentElement->Calculate(ARTIFICIAL_DIFFUSION_MATRIX, r_artificial_diffusion_matrix, rCurrentProcessInfo);

        LocalSystemMatrixType identity = identity_matrix<double>(rLHS_Contribution.size1());
        double coff = mStabilizationSourceCoefficient*mOverallDiffusionCoefficient*gauss_integration_weight;

        r_artificial_diffusion_matrix += coff*identity;

        noalias(rLHS_Contribution) -= r_artificial_diffusion_matrix;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{
        std::vector< LocalSystemMatrixType > mArtificialDiffusionMatrix;

        AdjointArtificialDiffusion mAdjointArtificialDiffusion;
        double mOverallDiffusionCoefficient;
        double mStabilizationSourceCoefficient;
        bool mIsMatrixEnergiesCalculated;
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedAdjointBossakScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED defined */
