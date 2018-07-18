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
#ifdef EIGEN_ROOT
    #include <Eigen/SVD>
    #include <Eigen/Eigenvalues>
#endif

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
#ifdef EIGEN_ROOT
        singularValuePressureCoupled,
        fullMatrixEigen,
#endif
        dynamicFullVMSSteadyMatrix,
        dynamicFullMatrix,
        deltaEnergyFullMatrixFullTimeStepped,
        deltaEnergyFullMatrixPartialTimeStepped
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
            "calculation_step" : 0
        })");

        rDiffusionParameters["stabilization_settings"].ValidateAndAssignDefaults(default_params);

        mEpsilon = rDiffusionParameters["stabilization_source_coefficient"].GetDouble();
        mTimeStep =  rDiffusionParameters["stabilization_settings"]["calculation_step"].GetInt();

        std::string method_name = rDiffusionParameters["stabilization_settings"]["method"].GetString();

        if (method_name=="dynamic_full_vms_steady_element_matrix")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::dynamicFullVMSSteadyMatrix;
        else if (method_name=="dynamic_full_element_matrix")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::dynamicFullMatrix;
        else if (method_name=="delta_energy_full_matrix_full_time_stepped")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::deltaEnergyFullMatrixFullTimeStepped;
        else if (method_name=="delta_energy_full_matrix_partial_time_stepped")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::deltaEnergyFullMatrixPartialTimeStepped;
#ifdef EIGEN_ROOT
        else if (method_name=="singular_value_pressure_coupled") {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::singularValuePressureCoupled;
            rDiffusionParameters["stabilization_source_coefficient"].SetDouble(0.0);
        }
        else if (method_name=="full_matrix_eigen_method")
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::fullMatrixEigen;
#endif
        else
        {
            if (method_name=="singular_value_pressure_coupled")
                KRATOS_ERROR<<"\"singular_value_pressure_coupled\" stabilization method is not compiled. Please compile Kratos with Eigen libraries.\n"<<rDiffusionParameters.PrettyPrintJsonString();
            else
                KRATOS_ERROR<<"stabilization method is not supported.\n"<<rDiffusionParameters.PrettyPrintJsonString();
        }

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
            case ArtificialDiffusionMethods::deltaEnergyFullMatrixFullTimeStepped:
                return CalculateArtificialDiffusionFullMatrixFullTimeStepped(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
            case ArtificialDiffusionMethods::deltaEnergyFullMatrixPartialTimeStepped:
                return CalculateArtificialDiffusionFullMatrixPartialTimeStepped(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
            case ArtificialDiffusionMethods::fullMatrixEigen:
                return CalculateArtificialDiffusionFullMatrixEigen(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
#ifdef EIGEN_ROOT
            case ArtificialDiffusionMethods::singularValuePressureCoupled:
                return CalculateArtificialDiffusionSVMethodPressureCoupled(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
#endif
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

#ifdef EIGEN_ROOT
    double CalculateVelocityDivergence(Element::Pointer pCurrentElement, ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType velocity_gradient;

        pCurrentElement->Calculate(VMS_VELOCITY_GRADIENT_TENSOR, velocity_gradient, rCurrentProcessInfo);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<velocity_gradient.size1(); i++)
            velocity_divergence += velocity_gradient(i,i);

        return velocity_divergence;
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 4, 4>  CalculateSVMethodCharacteristicMatrix(
                            Element::Pointer pCurrentElement,
                            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const Geometry< Node<3> >& r_geometry = pCurrentElement->GetGeometry();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();

        MatrixType velocity_gradient;

        pCurrentElement->Calculate(VMS_VELOCITY_GRADIENT_TENSOR, velocity_gradient, rCurrentProcessInfo);

        double velocity_divergence = 0.0;
        for (IndexType i=0;i<domain_size; i++)
            velocity_divergence += velocity_gradient(i,i);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 4, 4>  characteristic_matrix;
        characteristic_matrix.resize( domain_size+1, domain_size+1);
        for (IndexType i=0; i < domain_size; i++)
            characteristic_matrix(i,i) = 0.5 * velocity_divergence - velocity_gradient(i,i);
        for (IndexType i=0; i < domain_size; i++)
            for (IndexType j=i+1; j < domain_size; j++)
            {
                characteristic_matrix(i,j) =  velocity_gradient(i,j);
                characteristic_matrix(j,i) =  velocity_gradient(j,i);
            }

        return characteristic_matrix;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionSVMethodPressureCoupled(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const Geometry< Node<3> >& r_geometry = pCurrentElement->GetGeometry();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();

        auto  characteristic_matrix = CalculateSVMethodCharacteristicMatrix(pCurrentElement, rCurrentProcessInfo);

        for (IndexType i=0; i < domain_size + 1; i++)
        {
            characteristic_matrix(domain_size, i) = 0.0;
            characteristic_matrix(i, domain_size) = 0.0;
        }

        double velocity_divergence = CalculateVelocityDivergence(pCurrentElement, rCurrentProcessInfo);
        characteristic_matrix(domain_size, domain_size) = 0.5 * velocity_divergence;

        Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 4, 4>> svd(
                characteristic_matrix,
                Eigen::ComputeThinU | Eigen::ComputeThinV
                );

        const auto& S = svd.singularValues();

        double volume = 1.0;
        if (domain_size == 2)
            volume = pCurrentElement->GetGeometry().Area();
        else if (domain_size == 3)
            volume = pCurrentElement->GetGeometry().Volume();

        double artificial_diffusion = S[0];

        return artificial_diffusion;

        KRATOS_CATCH("");
    }

    double CalculateMaxEigenValue( const MatrixType& rMatrix)
    {
        const unsigned int matrix_size = rMatrix.size1();

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  symmetric_matrix;
        symmetric_matrix.resize(matrix_size, matrix_size);

        for (unsigned int i = 0; i < matrix_size; i++)
            for (unsigned int j = i; j < matrix_size; j++)
            {
                double value = 0.5*rMatrix(i,j) + 0.5*rMatrix(j,i);
                symmetric_matrix(i,j) = value;
                symmetric_matrix(j,i) = value;
            }

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es(symmetric_matrix, false);
        const auto& values = es.eigenvalues().real();

        return values[values.size()-1];

    }

    double CalculateArtificialDiffusionFullMatrixEigen(
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

        double previous_max_eigen_value = CalculateMaxEigenValue(rLHS_Contribution);
        double previous_beta = 0.0;

        // KRATOS_WATCH(pCurrentElement->Id());

        if (previous_max_eigen_value < 0.0)
            return 0.0;

        double current_beta = 1.0;
        double current_max_eigen_value = 0.0;

        unsigned int max_iterations = 10;

        double eigen_tolerance = 1e-6;
        for (unsigned int itr=0; itr < max_iterations; ++itr)
        {
            current_max_eigen_value = CalculateMaxEigenValue(rLHS_Contribution - current_beta*numerical_diffusion_matrix);
            double m = (current_max_eigen_value - previous_max_eigen_value)/(current_beta-previous_beta);
            previous_beta = current_beta;

            // KRATOS_WATCH(std::abs(current_max_eigen_value - previous_max_eigen_value));

            if (std::abs(current_max_eigen_value - previous_max_eigen_value) < eigen_tolerance)
                break;

            previous_max_eigen_value = current_max_eigen_value;
            current_beta -= current_max_eigen_value/m;
            // std::cout<<"Itr: "<<itr<<", beta: "<<previous_beta<<", max_eigen: "<<previous_max_eigen_value<<std::endl;
        }

        pCurrentElement->SetValue(ADJOINT_ENERGY, current_max_eigen_value);

        if (previous_beta > 0.0)
            return previous_beta;
        else
            return 0.0;

        KRATOS_CATCH("");
    }
#endif

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

        // double volume = 1.0;
        // if (domain_size == 2)
        //     volume = pCurrentElement->GetGeometry().Area();
        // else if (domain_size == 3)
        //     volume = pCurrentElement->GetGeometry().Volume();

        // artificial_diffusion /= volume;

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

        // double volume = 1.0;
        // if (domain_size == 2)
        //     volume = pCurrentElement->GetGeometry().Area();
        // else if (domain_size == 3)
        //     volume = pCurrentElement->GetGeometry().Volume();

        // artificial_diffusion /= volume;

        return artificial_diffusion;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionFullMatrixFullTimeStepped(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        return 0.0;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionFullMatrixPartialTimeStepped(
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
        double const adjoint_energy_current = inner_prod(temp_1, adjoint_values_vector);

        noalias(temp_1) = prod(numerical_diffusion_matrix, adjoint_values_vector);
        double const diffusion_energy = inner_prod(temp_1,adjoint_values_vector);

        pCurrentElement->GetValuesVector(adjoint_values_vector, mTimeStep + 1);
        noalias(temp_1) = prod(rLHS_Contribution, adjoint_values_vector);
        double const adjoint_energy_old = inner_prod(temp_1, adjoint_values_vector);

        double delta_energy = adjoint_energy_current - adjoint_energy_old;

        double artificial_diffusion = 0.0;

        KRATOS_DEBUG_ERROR_IF(diffusion_energy < 0.0)<<" --- Diffusion energy cannot be negative or zero."<<std::endl;

        if (delta_energy > 0.0)
        {
            artificial_diffusion = delta_energy/diffusion_energy;
        }

        pCurrentElement->SetValue(ADJOINT_ENERGY, delta_energy);
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
