//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 KratosFluidDynamicsApplication/license.txt
//
//  Main author:    Suneth Warnakulasriya, https://github.com/sunethwarna
//


#if !defined(KRATOS_ADJOINT_ARTIFICIAL_DIFFUSION)
#define KRATOS_ADJOINT_ARTIFICIAL_DIFFUSION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_elements/vms_adjoint_element.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// An Artificial Diffusion Calculation class to stabilize adjoint transient sensitivities.
class AdjointArtificialDiffusion
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointArtificialDiffusion);

    typedef Element::MatrixType MatrixType;

    typedef Element::GeometryType GeometryType;

    typedef Element::VectorType VectorType;

    // TODO: To be removed once diffusion methods are validated. Only one diffusion method will be kept
    enum class ArtificialDiffusionMethods
    {
        dynamicFullVMSSteadyMatrix,
        dynamicFullMatrix
    };

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    void SetArtificialDiffusionParameters(Parameters& rDiffusionParameters)
    {
        KRATOS_TRY;

        Parameters default_params(R"(
        {
            "method"    : "PLEASE_SPECIFY_A_METHOD",
            "time_step" : 0
        })");

        rDiffusionParameters["stabilization_settings"].ValidateAndAssignDefaults(default_params);

        mEpsilon = rDiffusionParameters["stabilization_source_coefficient"].GetDouble();
        mTimeStep =  rDiffusionParameters["stabilization_settings"]["time_step"].GetInt();

        std::string method_name = rDiffusionParameters["stabilization_settings"]["method"].GetString();

        if (method_name=="dynamic_full_vms_steady_element_matrix")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::dynamicFullVMSSteadyMatrix;
        else if (method_name=="dynamic_full_element_matrix")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::dynamicFullMatrix;
        else
            KRATOS_ERROR<<"stabilization method only supports \"dynamic_full_vms_steady_element_matrix\" or \"dynamic_full_element_matrix\""<<rDiffusionParameters.PrettyPrintJsonString();

        KRATOS_CATCH("");

    }

    double CalculateArtificialDiffusion(Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        switch (mArtificialDiffusionMethod)
        {
            case ArtificialDiffusionMethods::dynamicFullVMSSteadyMatrix:
                return CalculateArtificialDiffusionDynamicFullVMSSteadyMatrix(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
            case ArtificialDiffusionMethods::dynamicFullMatrix:
                return CalculateArtificialDiffusionDynamicFullMatrix(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
            default:
                return 0.0;
        }

        KRATOS_CATCH("");
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    ArtificialDiffusionMethods mArtificialDiffusionMethod = ArtificialDiffusionMethods::dynamicFullMatrix;
    double mEpsilon = 1e-6;
    int mTimeStep = 1;

    ///@}
    ///@name Private Operators
    ///@{

    double CalculateArtificialDiffusionDynamicFullVMSSteadyMatrix(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
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

        Matrix vms_steady_term_primal_gradient;
        pCurrentElement->Calculate(VMS_STEADY_TERM_PRIMAL_GRADIENT_MATRIX,
                                   vms_steady_term_primal_gradient,
                                   rCurrentProcessInfo);

        MatrixType numerical_diffusion_matrix;
        pCurrentElement->SetValue(ARTIFICIAL_DIFFUSION, 1.0);
        pCurrentElement->Calculate(ARTIFICIAL_DIFFUSION_MATRIX, numerical_diffusion_matrix, rCurrentProcessInfo);

        Matrix identity = identity_matrix<double>(vms_steady_term_primal_gradient.size1());
        numerical_diffusion_matrix += mEpsilon*gauss_integration_weight*identity;

        Vector adjoint_values_vector;
        pCurrentElement->GetValuesVector(adjoint_values_vector, mTimeStep);

        Vector temp_1;
        temp_1.resize(adjoint_values_vector.size());

        noalias(temp_1) = prod(vms_steady_term_primal_gradient, adjoint_values_vector);
        double const adjoint_energy = inner_prod(temp_1, adjoint_values_vector);

        noalias(temp_1) = prod(numerical_diffusion_matrix, adjoint_values_vector);
        double const diffusion_energy = inner_prod(temp_1,adjoint_values_vector);

        double artificial_diffusion = 0.0;

        KRATOS_DEBUG_ERROR_IF(diffusion_energy < 0.0)<<" --- Diffusion energy cannot be negative or zero."<<std::endl;

        if (adjoint_energy > 0.0)
        {
            artificial_diffusion = adjoint_energy/diffusion_energy;
        }

        pCurrentElement->SetValue(ADJOINT_ENERGY, adjoint_energy);
        pCurrentElement->SetValue(DIFFUSION_ENERGY, diffusion_energy);

        return artificial_diffusion;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionDynamicFullMatrix(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
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

        MatrixType numerical_diffusion_matrix;
        pCurrentElement->SetValue(ARTIFICIAL_DIFFUSION, 1.0);
        pCurrentElement->Calculate(ARTIFICIAL_DIFFUSION_MATRIX, numerical_diffusion_matrix, rCurrentProcessInfo);

        Matrix identity = identity_matrix<double>(rLHS_Contribution.size1());
        numerical_diffusion_matrix += mEpsilon*gauss_integration_weight*identity;

        Vector adjoint_values_vector;
        pCurrentElement->GetValuesVector(adjoint_values_vector, mTimeStep);

        Vector temp_1;
        temp_1.resize(adjoint_values_vector.size());

        noalias(temp_1) = prod(rLHS_Contribution, adjoint_values_vector);
        double const adjoint_energy = inner_prod(temp_1, adjoint_values_vector);

        noalias(temp_1) = prod(numerical_diffusion_matrix, adjoint_values_vector);
        double const diffusion_energy = inner_prod(temp_1,adjoint_values_vector);

        double artificial_diffusion = 0.0;

        KRATOS_DEBUG_ERROR_IF(diffusion_energy < 0.0)<<" --- Diffusion energy cannot be negative or zero."<<std::endl;

        if (adjoint_energy > 0.0)
        {
            artificial_diffusion = adjoint_energy/diffusion_energy;
        }

        pCurrentElement->SetValue(ADJOINT_ENERGY, adjoint_energy);
        pCurrentElement->SetValue(DIFFUSION_ENERGY, diffusion_energy);

        return artificial_diffusion;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_NUMERICAL_DIFFUSION defined */
