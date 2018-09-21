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
#include "utilities/svd_utils.h"

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
        VMSSteadyTermMatrixEigen,
        EnergyGenerationRateMatrixEigen,
#endif
        dynamicFullVMSSteadyMatrix,
        dynamicFullMatrix,
        energyGenerationRateWithAdjointVariables,
        energyGenerationRateMatrix
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

        KRATOS_INFO("Adjoint Stabilization")<<"--- Using stabilization source coefficient of "<<mEpsilon<<"."<<std::endl;
        KRATOS_INFO("Adjoint Stabilization")<<"--- Using time step "<<mTimeStep<<"."<<std::endl;

        std::string method_name = rDiffusionParameters["stabilization_settings"]["method"].GetString();

        if (method_name=="dynamic_full_vms_steady_element_matrix")
        {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::dynamicFullVMSSteadyMatrix;
            KRATOS_INFO("Adjoint Stabilization")<<"--- Using VMS Steady Element matrix for stabilization."<<std::endl;
        }
        else if (method_name=="energy_generation_rate_matrix")
        {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::energyGenerationRateMatrix;
            KRATOS_INFO("Adjoint Stabilization")<<"--- Using energy generation rate matrix for stabilization."<<std::endl;
        }
        else if (method_name=="energy_generation_rate_with_adjoint_variables")
        {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::energyGenerationRateWithAdjointVariables;
            KRATOS_INFO("Adjoint Stabilization")<<"--- Using energy generation rate matrix with adjoint variables for stabilization."<<std::endl;
        }
#ifdef EIGEN_ROOT
        else if (method_name=="singular_value_pressure_coupled")
        {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::singularValuePressureCoupled;
            rDiffusionParameters["stabilization_source_coefficient"].SetDouble(0.0);
            KRATOS_INFO("Adjoint Stabilization")<<"--- Using singular value method with pressure coupling for stabilization."<<std::endl;
        }
        else if (method_name=="vms_steady_term_matrix_eigen")
        {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::VMSSteadyTermMatrixEigen;
            KRATOS_INFO("Adjoint Stabilization")<<"--- Using VMS Steady Element matrix eigen values for stabilization."<<std::endl;
        }
        else if (method_name=="energy_generation_rate_matrix_eigen")
        {
            mArtificialDiffusionMethod = ArtificialDiffusionMethods::EnergyGenerationRateMatrixEigen;
            KRATOS_INFO("Adjoint Stabilization")<<"--- Using Energy generation rate matrix eigen values for stabilization."<<std::endl;
        }
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
            case ArtificialDiffusionMethods::energyGenerationRateWithAdjointVariables:
                return CalculateArtificialDiffusionEnergyGenerationRateWithAdjointVariables(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
            case ArtificialDiffusionMethods::energyGenerationRateMatrix:
                return CalculateArtificialDiffusionEnergyGenerationRateMatrix(
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
            case ArtificialDiffusionMethods::VMSSteadyTermMatrixEigen:
                return CalculateArtificialDiffusionVMSSteadyTermMatrixEigen(
                                        pCurrentElement,
                                        rLHS_Contribution,
                                        rRHS_Contribution,
                                        rCurrentProcessInfo);
            case ArtificialDiffusionMethods::EnergyGenerationRateMatrixEigen:
                return CalculateArtificialDiffusionEnergyGenerationRateMatrixEigen(
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

    double CalculateArtificialDiffusionBasedOnMatrixMaxEigenValue(
                                        Element::Pointer pCurrentElement,
                                        Matrix const& custom_matrix)
    {
        KRATOS_TRY;

        const Geometry< Node<3> >& r_geometry = pCurrentElement->GetGeometry();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();

        double gauss_integration_weight = 0.0;

        if (domain_size == 2)
            gauss_integration_weight = r_geometry.Area();
        else if (domain_size == 3)
            gauss_integration_weight = r_geometry.Volume();

        double custom_matrix_max_eigen_value = CalculateMaxEigenValue(custom_matrix);

        pCurrentElement->SetValue(STABILIZATION_ANALYSIS_MATRIX_CUSTOM, custom_matrix);
        pCurrentElement->SetValue(ADJOINT_ENERGY, custom_matrix_max_eigen_value);

        double numerical_diffusion = 0.0;
        if (custom_matrix_max_eigen_value > 0)
            numerical_diffusion = custom_matrix_max_eigen_value/(mEpsilon*gauss_integration_weight);

        return numerical_diffusion;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionEnergyGenerationRateMatrixEigen(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        Matrix vms_steady_term_primal_gradient;
        pCurrentElement->Calculate(VMS_STEADY_TERM_PRIMAL_GRADIENT_MATRIX,
                                   vms_steady_term_primal_gradient,
                                   rCurrentProcessInfo);

        return CalculateArtificialDiffusionBasedOnMatrixMaxEigenValue(pCurrentElement, vms_steady_term_primal_gradient);

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionVMSSteadyTermMatrixEigen(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        Matrix adjoint_energy_generation_matrix;
        pCurrentElement->Calculate(VMS_ADJOINT_ENERGY_GENERATION_RATE_MATRIX,
                                   adjoint_energy_generation_matrix,
                                   rCurrentProcessInfo);

        return CalculateArtificialDiffusionBasedOnMatrixMaxEigenValue(pCurrentElement, adjoint_energy_generation_matrix);

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

        KRATOS_ERROR_IF(diffusion_energy < 0.0)<<" --- Diffusion energy cannot be negative or zero."<<std::endl;

        if (std::abs(adjoint_energy) > 1e-3)
        {
            artificial_diffusion = std::abs(adjoint_energy/diffusion_energy);
        }

        pCurrentElement->SetValue(ADJOINT_ENERGY, adjoint_energy);
        pCurrentElement->SetValue(DIFFUSION_ENERGY, diffusion_energy);

        return artificial_diffusion;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionEnergyGenerationRateWithAdjointVariables(
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

        Matrix adjoint_energy_generation_matrix;
        pCurrentElement->Calculate(VMS_ADJOINT_ENERGY_GENERATION_RATE_MATRIX, adjoint_energy_generation_matrix, rCurrentProcessInfo);
        pCurrentElement->SetValue(STABILIZATION_ANALYSIS_MATRIX_CUSTOM, adjoint_energy_generation_matrix);

        noalias(temp_1) = prod(adjoint_energy_generation_matrix, adjoint_values_vector);
        double const adjoint_energy_rate = inner_prod(temp_1, adjoint_values_vector);

        noalias(temp_1) = prod(numerical_diffusion_matrix, adjoint_values_vector);
        double const diffusion_energy = inner_prod(temp_1,adjoint_values_vector);

        double artificial_diffusion = 0.0;

        KRATOS_ERROR_IF(diffusion_energy < 0.0)<<" --- Diffusion energy cannot be negative or zero."<<std::endl;

        if (adjoint_energy_rate > 0.0)
            artificial_diffusion = adjoint_energy_rate/diffusion_energy;

        pCurrentElement->SetValue(ADJOINT_ENERGY, adjoint_energy_rate);
        pCurrentElement->SetValue(DIFFUSION_ENERGY, diffusion_energy);

        return artificial_diffusion;

        KRATOS_CATCH("");
    }

    double CalculateArtificialDiffusionEnergyGenerationRateMatrix(
                                        Element::Pointer pCurrentElement,
                                        MatrixType& rLHS_Contribution,
                                        VectorType& rRHS_Contribution,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        MatrixType svd_u;
        MatrixType svd_v;
        MatrixType svd_s;

        unsigned int m_size = rLHS_Contribution.size1();

        svd_u.resize( m_size, m_size, false);
        svd_v.resize( m_size, m_size, false);
        svd_s.resize( m_size, m_size, false);

        svd_u.clear();
        svd_v.clear();
        svd_s.clear();

        Matrix adjoint_energy_generation_matrix;
        pCurrentElement->Calculate(VMS_ADJOINT_ENERGY_GENERATION_RATE_MATRIX, adjoint_energy_generation_matrix, rCurrentProcessInfo);

        SVDUtils<double>::SingularValueDecomposition(adjoint_energy_generation_matrix, svd_u, svd_s, svd_v);

        return svd_s(0,0);

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
