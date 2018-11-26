//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_REYNOLDS_STRESS_TENSOR_H_INCLUDED)
#define KRATOS_REYNOLDS_STRESS_TENSOR_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <string>

// Project includes
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
class ReynoldsStressTensor
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ReynoldsStressTensor);

    constexpr static unsigned int TBlockSize = TDim + 1;

    constexpr static unsigned int TFluidLocalSize = TBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    constexpr static double INV_TDIM = 1.0 / TDim;

    ReynoldsStressTensor(BoundedMatrix<double, TNumNodes, TDim>& rVelocity,
                         BoundedMatrix<double, TNumNodes, TDim>& rIntegrationMatrix,
                         Vector& rCoefficients,
                         Geometry<Node<3>>& rGeometry,
                         double TurbulentKineticEnergy,
                         double TurbulentKinematicViscosity,
                         double Density,
                         double GaussPointWeight)
        : mVelocity(rVelocity), mIntegrationMatrix(rIntegrationMatrix), mCoefficients(rCoefficients), mGeometry(rGeometry)
    {
        mTurbulentKineticEnergy = TurbulentKineticEnergy;
        mTurbulentKinematicViscosity = TurbulentKinematicViscosity;
        mDensity = Density;
        mGaussPointWeight = GaussPointWeight;
    }

    void test()
    {
        ShapeParameter Deriv;
        Deriv.NodeIndex = 2;
        Deriv.Direction = 0;
        const double delta = 1e-12;

        BoundedMatrix<double, TDim, TDim> matrix_0;
        BoundedMatrix<double, TDim, TDim> matrix_1;
        BoundedMatrix<double, TDim, TDim> matrix_fd;
        BoundedMatrix<double, TDim, TDim> analytical_matrix_output;

        Vector coffs;
        coffs.resize(9);
        coffs[0] = 3.0;
        coffs[1] = 5.0;
        coffs[2] = 7.0;
        coffs[4] = 11.0;
        coffs[5] = 13.0;
        coffs[6] = 17.0;
        coffs[7] = 19.0;
        coffs[8] = 23.0;

        mCoefficients = coffs;

        std::cout<<"1.1"<<std::endl;
        CalculateReynoldsStressTensor(matrix_0);
        std::cout<<"1.2"<<std::endl;        
        CalculateReynoldsStressTensorPrimalDerivative(Deriv, analytical_matrix_output);
        std::cout<<"1.3"<<std::endl;
        mVelocity(Deriv.NodeIndex, Deriv.Direction) += delta;
        std::cout<<"1.4"<<std::endl;
        CalculateReynoldsStressTensor(matrix_1);
        std::cout<<"1.5"<<std::endl;
        mVelocity(Deriv.NodeIndex, Deriv.Direction) -= delta;
        std::cout<<"1.6"<<std::endl;
        CalculateFDMatrixSensitivity(matrix_0, matrix_1, delta, matrix_fd);
        std::cout<<"1.7"<<std::endl;
        AssertMatrices(analytical_matrix_output, matrix_fd);
        std::cout<<"1.8"<<std::endl;

        array_1d<double, TNumNodes> N;
        std::cout<<"1.9"<<std::endl;
        CalculateReynoldsStressTensorShapeDerivative(Deriv, analytical_matrix_output);
        std::cout<<"1.10"<<std::endl;
        const double edge_length = mGeometry.MinEdgeLength()*1e-8;
        std::cout<<"1.11"<<std::endl;
        mGeometry[Deriv.NodeIndex].Coordinates()[Deriv.Direction] += edge_length;
        std::cout<<"1.12"<<std::endl;
        GeometryUtils::CalculateGeometryData(mGeometry, mIntegrationMatrix, N, mGaussPointWeight);
        std::cout<<"1.13"<<std::endl;
        CalculateReynoldsStressTensor(matrix_1);
        std::cout<<"1.14"<<std::endl;
        mGeometry[Deriv.NodeIndex].Coordinates()[Deriv.Direction] -= edge_length;
        std::cout<<"1.15"<<std::endl;
        CalculateFDMatrixSensitivity(matrix_0, matrix_1, edge_length, matrix_fd);
        std::cout<<"1.16"<<std::endl;
        AssertMatrices(analytical_matrix_output, matrix_fd);
        std::cout<<"1.17"<<std::endl;

        // test sensitivity of symmetric gradient of velocity gradient
        // CalculateSymmetricVelocityGradientShapeDerivative(Deriv, analytical_matrix_output);
        // CalculateSymmetricVelocityGradient(velocity_0, matrix_0);
        // auto& rGeom = mrCurrentElement.GetGeometry();
        // double edge_length = rGeom.MinEdgeLength()*1e-8;
        // rGeom[Deriv.NodeIndex].Coordinates()[Deriv.Direction] += edge_length;
        // GeometryUtils::CalculateGeometryData(mrCurrentElement.GetGeometry(), mIntegrationMatrix, mN, mGaussPointWeight);
        // CalculateSymmetricVelocityGradient(velocity_0, matrix_1);
        // CalculateFDMatrixSensitivity(matrix_0, matrix_1, edge_length, matrix_fd);
        // AssertMatrices(analytical_matrix_output, matrix_fd);


        // CalculateAntiSymmetricVelocityGradientShapeDerivative(Deriv, analytical_matrix_output);
        // CalculateAntiSymmetricVelocityGradient(velocity_0, matrix_0);
        // auto& rGeom = mrCurrentElement.GetGeometry();
        // double edge_length = rGeom.MinEdgeLength()*1e-8;
        // rGeom[Deriv.NodeIndex].Coordinates()[Deriv.Direction] += edge_length;
        // GeometryUtils::CalculateGeometryData(mrCurrentElement.GetGeometry(), mIntegrationMatrix, mN, mGaussPointWeight);
        // CalculateAntiSymmetricVelocityGradient(velocity_0, matrix_1);
        // CalculateFDMatrixSensitivity(matrix_0, matrix_1, edge_length, matrix_fd);
        // AssertMatrices(analytical_matrix_output, matrix_fd);

        // CalculateAntiSymmetricVelocityGradientPrimalDerivative(Deriv,
        // analytical_matrix_output); CalculateAntiSymmetricVelocityGradient(velocity_0,
        // matrix_0); CalculateAntiSymmetricVelocityGradient(velocity_1,
        // matrix_1); CalculateFDMatrixSensitivity(matrix_0, matrix_1, delta,
        // matrix_fd); AssertMatrices(analytical_matrix_output, matrix_fd);

        // double norm_0_val = norm_frobenius(matrix_0);
        // norm_0_val *= norm_0_val;
        // double norm_1_val = norm_frobenius(matrix_1);
        // norm_1_val *= norm_1_val;

        // double fd_sensitivity = (norm_1_val - norm_0_val)/delta;
        // double analytical_sensitivity;
        // CalculateFrobeniusNormSquareDerivative(matrix_0, analytical_matrix_output, analytical_sensitivity);
        // if (std::abs((analytical_sensitivity-fd_sensitivity)/fd_sensitivity) > 1e-3)
        //     std::cout<<std::scientific<<fd_sensitivity<<", "<<analytical_sensitivity<<std::endl;
    }

    void CalculateFDMatrixSensitivity(const BoundedMatrix<double, TDim, TDim> matrix_0,
                                      const BoundedMatrix<double, TDim, TDim> matrix_1,
                                      const double epsilon,
                                      BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput = matrix_1 - matrix_0;
        rOutput *= (1.0 / epsilon);
    }

    void AssertMatrices(const BoundedMatrix<double, TDim, TDim> matrix_0,
                        BoundedMatrix<double, TDim, TDim> matrix_1)
    {
        KRATOS_WATCH(matrix_0);
        KRATOS_WATCH(matrix_1);
        const double tol = 2e-3;
        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
            {
                if (matrix_1(i, j) == 0)
                {
                    if (std::abs(matrix_0(i, j)) > 1e-6)
                    {
                        std::cout << "Test Failed." << std::endl;
                        std::cout << "Matrix 0:" << matrix_0 << std::endl;
                        std::cout << "Matrix 1:" << matrix_1 << std::endl;

                        std::cout << std::scientific << matrix_0(i, j)
                                  << "!=" << matrix_1(i, j) << std::endl;

                        std::exit(-1);
                    }
                }
                else
                {
                    const double relative_tolerance =
                        std::abs((matrix_0(i, j) - matrix_1(i, j)) / matrix_1(i, j));
                    if (relative_tolerance > tol)
                    {
                        std::cout << "Test Failed." << std::endl;
                        std::cout << "Matrix 0:" << matrix_0 << std::endl;
                        std::cout << "Matrix 1:" << matrix_1 << std::endl;

                        std::cout << std::scientific << matrix_0(i, j)
                                  << "!=" << matrix_1(i, j)
                                  << "[ difference = " << relative_tolerance
                                  << "]" << std::endl;
                        std::exit(-1);
                    }
                }
            }
    }

    void AddReynoldsStressTensorVelocityContributionLHS(Matrix& rLeftHandSideMatrix)
    {
        KRATOS_TRY

        ShapeParameter Deriv;

        BoundedMatrix<double, TDim, TDim> SymmetricVelGrad;
        CalculateSymmetricVelocityGradient(mVelocity, SymmetricVelGrad);

        // Calculating left hand side matrix - only the linear part is added to the LHS
        for (unsigned int b = 0; b < TNumNodes; ++b)
        {
            for (unsigned int j = 0; j < TDim; ++j)
            {

                Deriv.NodeIndex = b;
                Deriv.Direction = j;

                BoundedMatrix<double, TDim, TDim> SymmetricVelGradDerivative;
                CalculateSymmetricVelocityGradientPrimalDerivative(Deriv, SymmetricVelGradDerivative);

                const double SymmetricVelGradDerivTrace = CalculateMatrixTrace(SymmetricVelGradDerivative);

                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int i = 0; i < TDim; ++i)
                    {

                        double valij = 0;
                        for (unsigned int l = 0; l < TDim; ++l)
                            valij += mTurbulentKinematicViscosity * SymmetricVelGradDerivative(l,i) * mIntegrationMatrix(a,l);

                        valij -= mTurbulentKinematicViscosity * INV_TDIM * mIntegrationMatrix(a,i) * SymmetricVelGradDerivTrace;
                        
                        rLeftHandSideMatrix(a*TBlockSize+i, b*TBlockSize+j) += mDensity * valij * mGaussPointWeight;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void AddReynoldsStressTensorVelocityContributionRHS(Vector& rRightHandSideVector)
    {
        KRATOS_TRY

        BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
        reynolds_stress_tensor.clear();
        AddHydrostaticReynoldsStressTensor(reynolds_stress_tensor);
        AddDeviatoricReynoldsStressTensorNonLinearPart(reynolds_stress_tensor);

        for (unsigned int a = 0; a < TNumNodes; ++a)
        {
            for (unsigned int i = 0; i < TDim; ++i)
            {
                for (unsigned int l = 0; l < TDim; ++l)
                    rRightHandSideVector(a * TBlockSize + i) -= mDensity * mIntegrationMatrix(a,l) * reynolds_stress_tensor(l,i) * mGaussPointWeight;
            }
        }
        KRATOS_CATCH("")
    }


    void AddReynoldsStressTensorPrimalDerivativeContributionLHS(Matrix& rLeftHandSideMatrix,
                                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ShapeParameter Deriv;
        BoundedMatrix<double, TDim, TDim> ReynoldsStressTensorDerivative;

        for (unsigned int c = 0; c < TNumNodes; ++c)
        {
            for (unsigned int k = 0; k < TDim; ++k)
            {
                this->CalculateReynoldsStressTensorPrimalDerivative(Deriv, ReynoldsStressTensorDerivative);
                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int i = 0; i < TDim; ++i)
                    {
                        double valik = 0.0;
                        for (unsigned int l = 0; l < TDim; ++l)
                        {
                            valik += mIntegrationMatrix(a,l) * ReynoldsStressTensorDerivative(l, i);
                        }

                        rLeftHandSideMatrix(a*TBlockSize+i, c*TBlockSize + k) += mDensity * valik * mGaussPointWeight;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void AddReynoldsStressTensorShapeDerivativeContribution(Matrix& rShapeDerivativesMatrix,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
        ShapeParameter Deriv;
        BoundedMatrix<double, TDim, TDim> ReynoldsStressTensorValue;
        BoundedMatrix<double, TDim, TDim> ReynoldsStressTensorDerivative;
        Matrix DN_DX_Deriv;
        double DetJ0_deriv;

        this->CalculateReynoldsStressTensor(ReynoldsStressTensorValue);

        for (unsigned int r = 0; r < TNumNodes; ++r)
        {
            for (unsigned int m = 0; m < TDim; ++m)
            {
                Deriv.NodeIndex = r;
                Deriv.Direction = m;
                this->CalculateReynoldsStressTensorShapeDerivative(Deriv, ReynoldsStressTensorDerivative);
                this->CalculateShapeFunctionDerivativeShapeDerivatives(0, Deriv, DetJ0_deriv, DN_DX_Deriv);

                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int i = 0; i < TDim; ++i)
                    {
                        double valim = 0.0;
                        for (unsigned int l = 0; l < TDim; ++l)
                        {
                            valim += mIntegrationMatrix(a,l) * ReynoldsStressTensorValue(l,i) * DetJ0_deriv;
                            valim += mIntegrationMatrix(a,l) * ReynoldsStressTensorDerivative(l,i) * mGaussPointWeight;
                            valim += DN_DX_Deriv(a,l) * ReynoldsStressTensorValue(l,i) * mGaussPointWeight;
                        }

                        rShapeDerivativesMatrix(r * TCoordLocalSize + m, a * TBlockSize + i) += mDensity * valim;
                    }
                }
            }
        }
    }

    void AddReynoldsStressTensorConditionVelocityContribution()
    {

    }

    void AddReynoldsStressTensorPrimalDerivativeConditionContribution()
    {

    }

    void AddReynoldsStressTensorShapeDerivativeConditionContribution()
    {

    }

    void CalculateReynoldsStressTensor(BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        AddDeviatoricReynoldsStressTensorLinearPart(mVelocity, rOutput);
        AddDeviatoricReynoldsStressTensorNonLinearPart(mVelocity, rOutput);
        AddHydrostaticReynoldsStressTensor(rOutput);
    }

    void CalculateReynoldsStressTensorPrimalDerivative(const ShapeParameter& Deriv,
                                                       BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        CalculateReynoldsStressTensorPrimalDerivative(Deriv, mVelocity, rOutput);
    }

    void CalculateReynoldsStressTensorShapeDerivative(const ShapeParameter& Deriv,
                                                      BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        CalculateReynoldsStressTensorShapeDerivative(Deriv, mVelocity, rOutput);
    }

private:
    ///@name Member Variables
    ///@{
    BoundedMatrix<double, TNumNodes, TDim>& mVelocity;
    BoundedMatrix<double, TNumNodes, TDim>& mIntegrationMatrix;
    Vector& mCoefficients;
    Geometry<Node<3>>& mGeometry;
    double mDensity;
    double mGaussPointWeight;
    double mTurbulentKineticEnergy;
    double mTurbulentKinematicViscosity;

    double CalculateMatrixElementProductSum(BoundedMatrix<double, TDim, TDim> matrix_1, BoundedMatrix<double, TDim, TDim> matrix_2)
    {
        double result = 0.0;
        for (unsigned int i = 0; i < TDim; ++i)
            for (unsigned int j = 0; j < TDim; ++j)
                result += matrix_1(i,j) * matrix_2 (i,j);

        return result;
    }

    double CalculateMatrixTrace(const Matrix& rMatrix)
    {
        double result = 0.0;
        for (unsigned int i = 0; i < rMatrix.size1(); i++)
            result += rMatrix(i, i);
        return result;
    }

    void CalculateTwoMatrixProductDerivative(
        const BoundedMatrix<double, TDim, TDim>& rMatrix1,
        const BoundedMatrix<double, TDim, TDim>& rMatrix1Derivative,
        const BoundedMatrix<double, TDim, TDim>& rMatrix2,
        const BoundedMatrix<double, TDim, TDim>& rMatrix2Derivative,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        rOutput = prod(rMatrix1, rMatrix2Derivative) + prod(rMatrix1Derivative, rMatrix2);
    }

    void CalculateThreeMatrixProductDerivative(
        const BoundedMatrix<double, TDim, TDim>& rMatrix1,
        const BoundedMatrix<double, TDim, TDim>& rMatrix1Derivative,
        const BoundedMatrix<double, TDim, TDim>& rMatrix2,
        const BoundedMatrix<double, TDim, TDim>& rMatrix2Derivative,
        const BoundedMatrix<double, TDim, TDim>& rMatrix3,
        const BoundedMatrix<double, TDim, TDim>& rMatrix3Derivative,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        BoundedMatrix<double, TDim, TDim> temp;
        temp = prod(rMatrix2, rMatrix3);
        rOutput += prod(rMatrix1Derivative, temp);
        temp = prod(rMatrix2Derivative, rMatrix3);
        rOutput += prod(rMatrix1, temp);
        temp = prod(rMatrix2, rMatrix3Derivative);
        rOutput += prod(rMatrix1, temp);
    }

    void CalculateSymmetricVelocityGradientPrimalDerivative(
        const ShapeParameter Deriv, BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        const std::size_t c = Deriv.NodeIndex;
        const std::size_t k = Deriv.Direction;

        for (unsigned int i = 0; i < TDim; ++i)
            for (unsigned int j = 0; j < TDim; ++j)
            {
                if (i == k)
                    rOutput(i, j) += mIntegrationMatrix(c, j);
                if (j == k)
                    rOutput(i, j) += mIntegrationMatrix(c, i);
            }

        rOutput *= 0.5;

        KRATOS_CATCH("")
    }

    void CalculateShapeFunctionDerivativeShapeDerivatives(const unsigned int IntegrationPointIndex,
                                                          const ShapeParameter Deriv,
                                                          double& rDetJ0_deriv,
                                                          Matrix& rDN_DX0_deriv)
    {
        KRATOS_TRY

        const auto& r_geom = mGeometry;

        const auto integration_method =
            IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);

        Matrix J0(TDim, TDim);
        Matrix DN_DX0_deriv;
        const auto& integration_points = r_geom.IntegrationPoints(integration_method);

        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[IntegrationPointIndex], J0);
        const Matrix& rDN_De =
            r_geom.ShapeFunctionsLocalGradients(integration_method)[IntegrationPointIndex];

        GeometricalSensitivityUtility geometrical_sensitivity(J0, rDN_De);
        geometrical_sensitivity.CalculateSensitivity(Deriv, rDetJ0_deriv, rDN_DX0_deriv);

        KRATOS_CATCH("")
    }

    void CalculateSymmetricVelocityGradientShapeDerivative(const ShapeParameter Deriv, BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        Matrix rDN_DX0_deriv;
        double DetJ0_deriv;
        CalculateShapeFunctionDerivativeShapeDerivatives(0, Deriv, DetJ0_deriv, rDN_DX0_deriv);

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
            {
                for (unsigned int d = 0; d < TNumNodes; d++)
                {
                    rOutput(i,j) += mVelocity(d,i) * rDN_DX0_deriv(d,j);
                    rOutput(i,j) += mVelocity(d,j) * rDN_DX0_deriv(d,i);
                }
            }

        rOutput *= 0.5;

        KRATOS_CATCH("")
    }

    void CalculateAntiSymmetricVelocityGradientShapeDerivative(const ShapeParameter Deriv, BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        Matrix rDN_DX0_deriv;
        double DetJ0_deriv;
        CalculateShapeFunctionDerivativeShapeDerivatives(0, Deriv, DetJ0_deriv, rDN_DX0_deriv);

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
            {
                for (unsigned int d = 0; d < TNumNodes; d++)
                {
                    rOutput(i,j) += mVelocity(d,i) * rDN_DX0_deriv(d,j);
                    rOutput(i,j) -= mVelocity(d,j) * rDN_DX0_deriv(d,i);
                }
            }

        rOutput *= 0.5;

        KRATOS_CATCH("")
    }

    void CalculateAntiSymmetricVelocityGradientPrimalDerivative(
        const ShapeParameter Deriv, BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        const std::size_t c = Deriv.NodeIndex;
        const std::size_t k = Deriv.Direction;

        for (unsigned int i = 0; i < TDim; ++i)
            for (unsigned int j = 0; j < TDim; ++j)
            {
                if (i == k)
                    rOutput(i, j) += mIntegrationMatrix(c, j);
                if (j == k)
                    rOutput(i, j) -= mIntegrationMatrix(c, i);
            }

        rOutput *= 0.5;

        KRATOS_CATCH("")
    }

    double CalculateFrobeniusNormSquareDerivative(
        const BoundedMatrix<double, TDim, TDim>& rMatrix,
        const BoundedMatrix<double, TDim, TDim>& rMatrixDeriv)
    {
        double rOutput = 0.0;
        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                rOutput += rMatrix(i, j) * rMatrixDeriv(i, j);

        return rOutput * 2.0;
    }

    void CalculateSymmetricVelocityGradient(const BoundedMatrix<double, TNumNodes, TDim>& velocity,
                                            BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                for (unsigned int d = 0; d < TNumNodes; d++)
                {
                    rOutput(i, j) += mIntegrationMatrix(d, j) * velocity(d, i);
                    rOutput(i, j) += mIntegrationMatrix(d, i) * velocity(d, j);
                }

        rOutput *= 0.5;

        KRATOS_CATCH("")
    }

    void CalculateAntiSymmetricVelocityGradient(const BoundedMatrix<double, TNumNodes, TDim>& velocity,
                                                BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                for (unsigned int d = 0; d < TNumNodes; d++)
                {
                    rOutput(i, j) += mIntegrationMatrix(d, j) * velocity(d, i);
                    rOutput(i, j) -= mIntegrationMatrix(d, i) * velocity(d, j);
                }

        rOutput *= 0.5;

        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorNonLinearPart(const BoundedMatrix<double, TNumNodes, TDim>& rVelocity,
                                                        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);
        BoundedMatrix<double, TDim, TDim> temp;

        BoundedMatrix<double, TDim, TDim> SymmetricVelGrad;
        CalculateSymmetricVelocityGradient(rVelocity, SymmetricVelGrad);
        BoundedMatrix<double, TDim, TDim> AntiSymmetricVelGrad;
        CalculateAntiSymmetricVelocityGradient(rVelocity, AntiSymmetricVelGrad);

        const double symmetric_frobenius_norm_square =
            std::pow(norm_frobenius(SymmetricVelGrad), 2);
        const double anti_symmetric_frobenius_norm_square =
            std::pow(norm_frobenius(AntiSymmetricVelGrad), 2);

        double coefficient = mCoefficients[8] * std::pow(mCoefficients[0], 2);
        rOutput -= coefficient * mCoefficients[1] *
                   (prod(SymmetricVelGrad, SymmetricVelGrad) -
                    INV_TDIM * symmetric_frobenius_norm_square * mIdentity);
        rOutput -= coefficient * mCoefficients[2] * prod(AntiSymmetricVelGrad, SymmetricVelGrad);
        rOutput += coefficient * mCoefficients[2] * prod(SymmetricVelGrad, AntiSymmetricVelGrad);

        rOutput -= coefficient * mCoefficients[3] *
                   (prod(AntiSymmetricVelGrad, trans(AntiSymmetricVelGrad)) -
                    INV_TDIM * anti_symmetric_frobenius_norm_square * mIdentity);

        coefficient = mCoefficients[8] * std::pow(mCoefficients[0], 3);
        temp = prod(SymmetricVelGrad, AntiSymmetricVelGrad);
        rOutput -= coefficient * mCoefficients[4] * prod(SymmetricVelGrad, temp);
        temp = prod(SymmetricVelGrad, SymmetricVelGrad);
        rOutput += coefficient * mCoefficients[4] * prod(AntiSymmetricVelGrad, temp);

        temp = prod(AntiSymmetricVelGrad, AntiSymmetricVelGrad);
        rOutput -= coefficient * mCoefficients[5] * prod(temp, SymmetricVelGrad);
        rOutput += coefficient * mCoefficients[5] * (2 * INV_TDIM) *
                   CalculateMatrixTrace(prod(temp, SymmetricVelGrad)) * mIdentity;
        temp = prod(SymmetricVelGrad, AntiSymmetricVelGrad);
        rOutput -= coefficient * mCoefficients[5] * prod(temp, AntiSymmetricVelGrad);

        rOutput -= coefficient *
                               (mCoefficients[6] * symmetric_frobenius_norm_square +
                                mCoefficients[7] * anti_symmetric_frobenius_norm_square) *
            (SymmetricVelGrad - INV_TDIM * CalculateMatrixTrace(SymmetricVelGrad) * mIdentity);

        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorLinearPart(const BoundedMatrix<double, TNumNodes, TDim>& rVelocity,
                                                     BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);

        BoundedMatrix<double, TDim, TDim> SymmetricVelGrad;
        CalculateSymmetricVelocityGradient(rVelocity, SymmetricVelGrad);

        rOutput += mTurbulentKinematicViscosity * (SymmetricVelGrad - INV_TDIM * CalculateMatrixTrace(SymmetricVelGrad) * mIdentity);

        KRATOS_CATCH("")        
    }

    void AddHydrostaticReynoldsStressTensor(BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);        

        rOutput -= 2 * INV_TDIM * mTurbulentKineticEnergy * mIdentity;

        KRATOS_CATCH("")
    }

    void CalculateReynoldsStressTensorPrimalDerivative(
        const ShapeParameter& Deriv,
        const BoundedMatrix<double, TNumNodes, TDim>& rVelocity,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);
        BoundedMatrix<double, TDim, TDim> temp;

        rOutput.clear();

        BoundedMatrix<double, TDim, TDim> SymmetricVelGrad;
        CalculateSymmetricVelocityGradient(rVelocity, SymmetricVelGrad);
        BoundedMatrix<double, TDim, TDim> AntiSymmetricVelGrad;
        CalculateAntiSymmetricVelocityGradient(rVelocity, AntiSymmetricVelGrad);
        BoundedMatrix<double, TDim, TDim> SymmetricVelGradDerivative;
        CalculateSymmetricVelocityGradientPrimalDerivative(Deriv, SymmetricVelGradDerivative);
        BoundedMatrix<double, TDim, TDim> AntiSymmetricVelGradDerivative;
        CalculateAntiSymmetricVelocityGradientPrimalDerivative(
            Deriv, AntiSymmetricVelGradDerivative);

        const double symmetric_frobenius_norm_square =
            std::pow(norm_frobenius(SymmetricVelGrad), 2);
        const double anti_symmetric_frobenius_norm_square =
            std::pow(norm_frobenius(AntiSymmetricVelGrad), 2);
        const double symmetric_frobenius_norm_square_derivative =
            CalculateFrobeniusNormSquareDerivative(
                SymmetricVelGrad, SymmetricVelGradDerivative);
        const double anti_symmetric_frobenius_norm_square_derivative =
            CalculateFrobeniusNormSquareDerivative(
                AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative);

        double coefficient = -mCoefficients[8] * mCoefficients[0];

        rOutput +=
            (coefficient - coefficient * std::pow(mCoefficients[0], 2) *
                               (mCoefficients[6] * symmetric_frobenius_norm_square +
                                mCoefficients[7] * anti_symmetric_frobenius_norm_square)) *
            (SymmetricVelGradDerivative -
             INV_TDIM * mIntegrationMatrix(Deriv.NodeIndex, Deriv.Direction) * mIdentity);
        temp = SymmetricVelGrad - INV_TDIM * CalculateMatrixTrace(SymmetricVelGrad) * mIdentity;
        rOutput -= coefficient * std::pow(mCoefficients[0], 2) * mCoefficients[6] *
                   symmetric_frobenius_norm_square_derivative * temp;
        rOutput -= coefficient * std::pow(mCoefficients[0], 2) * mCoefficients[7] *
                   anti_symmetric_frobenius_norm_square_derivative * temp;

        coefficient = mCoefficients[8] * std::pow(mCoefficients[0], 2);

        CalculateTwoMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, SymmetricVelGrad,
            SymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[1] *
                   (temp - INV_TDIM * symmetric_frobenius_norm_square_derivative * mIdentity);

        CalculateTwoMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            SymmetricVelGrad, SymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[2] * (temp);
        CalculateTwoMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, temp);
        rOutput -= coefficient * mCoefficients[2] * (temp);

        CalculateTwoMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            -AntiSymmetricVelGrad, -AntiSymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[3] *
                   (temp - INV_TDIM * anti_symmetric_frobenius_norm_square_derivative * mIdentity);

        coefficient = mCoefficients[8] * std::pow(mCoefficients[0], 3);
        CalculateThreeMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, SymmetricVelGrad,
            SymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[4] * (temp);
        CalculateThreeMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            SymmetricVelGrad, SymmetricVelGradDerivative, SymmetricVelGrad,
            SymmetricVelGradDerivative, temp);
        rOutput -= coefficient * mCoefficients[4] * (temp);

        CalculateThreeMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            SymmetricVelGrad, SymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[5] * (temp);
        rOutput -= coefficient * mCoefficients[5] * 2 * INV_TDIM *
                   CalculateMatrixTrace(temp) * mIdentity;
        CalculateThreeMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[5] * (temp);

        rOutput *= -1.0;

        KRATOS_CATCH("")
    }

    void CalculateReynoldsStressTensorShapeDerivative(
        const ShapeParameter& Deriv,
        const BoundedMatrix<double, TNumNodes, TDim>& rVelocity,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);
        BoundedMatrix<double, TDim, TDim> temp;

        rOutput.clear();

        BoundedMatrix<double, TDim, TDim> SymmetricVelGrad;
        CalculateSymmetricVelocityGradient(rVelocity, SymmetricVelGrad);
        BoundedMatrix<double, TDim, TDim> AntiSymmetricVelGrad;
        CalculateAntiSymmetricVelocityGradient(rVelocity, AntiSymmetricVelGrad);
        BoundedMatrix<double, TDim, TDim> SymmetricVelGradDerivative;
        CalculateSymmetricVelocityGradientShapeDerivative(Deriv, SymmetricVelGradDerivative);
        BoundedMatrix<double, TDim, TDim> AntiSymmetricVelGradDerivative;
        CalculateAntiSymmetricVelocityGradientShapeDerivative(
            Deriv, AntiSymmetricVelGradDerivative);

        const double symmetric_frobenius_norm_square =
            std::pow(norm_frobenius(SymmetricVelGrad), 2);
        const double anti_symmetric_frobenius_norm_square =
            std::pow(norm_frobenius(AntiSymmetricVelGrad), 2);
        const double symmetric_frobenius_norm_square_derivative =
            CalculateFrobeniusNormSquareDerivative(
                SymmetricVelGrad, SymmetricVelGradDerivative);
        const double anti_symmetric_frobenius_norm_square_derivative =
            CalculateFrobeniusNormSquareDerivative(
                AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative);

        double coefficient = -mCoefficients[8] * mCoefficients[0];
        rOutput +=
            (coefficient - coefficient * std::pow(mCoefficients[0], 2) *
                               (mCoefficients[6] * symmetric_frobenius_norm_square +
                                mCoefficients[7] * anti_symmetric_frobenius_norm_square)) *
            (SymmetricVelGradDerivative -
             INV_TDIM * CalculateMatrixTrace(SymmetricVelGradDerivative) * mIdentity);
        temp = SymmetricVelGrad - INV_TDIM * CalculateMatrixTrace(SymmetricVelGrad) * mIdentity;
        rOutput -= coefficient * std::pow(mCoefficients[0], 2) * mCoefficients[6] *
                   symmetric_frobenius_norm_square_derivative * temp;
        rOutput -= coefficient * std::pow(mCoefficients[0], 2) * mCoefficients[7] *
                   anti_symmetric_frobenius_norm_square_derivative * temp;

        coefficient = mCoefficients[8] * std::pow(mCoefficients[0], 2);

        CalculateTwoMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, SymmetricVelGrad,
            SymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[1] *
                   (temp - INV_TDIM * symmetric_frobenius_norm_square_derivative * mIdentity);

        CalculateTwoMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            SymmetricVelGrad, SymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[2] * (temp);
        CalculateTwoMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, temp);
        rOutput -= coefficient * mCoefficients[2] * (temp);

        CalculateTwoMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            -AntiSymmetricVelGrad, -AntiSymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[3] *
                   (temp - INV_TDIM * anti_symmetric_frobenius_norm_square_derivative * mIdentity);

        coefficient = mCoefficients[8] * std::pow(mCoefficients[0], 3);
        CalculateThreeMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, SymmetricVelGrad,
            SymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[4] * (temp);
        CalculateThreeMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            SymmetricVelGrad, SymmetricVelGradDerivative, SymmetricVelGrad,
            SymmetricVelGradDerivative, temp);
        rOutput -= coefficient * mCoefficients[4] * (temp);

        CalculateThreeMatrixProductDerivative(
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            AntiSymmetricVelGrad, AntiSymmetricVelGradDerivative,
            SymmetricVelGrad, SymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[5] * (temp);
        rOutput -= coefficient * mCoefficients[5] * 2 * INV_TDIM *
                   CalculateMatrixTrace(temp) * mIdentity;
        CalculateThreeMatrixProductDerivative(
            SymmetricVelGrad, SymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, AntiSymmetricVelGrad,
            AntiSymmetricVelGradDerivative, temp);
        rOutput += coefficient * mCoefficients[5] * (temp);

        rOutput *= (-1.0);

        KRATOS_CATCH("")
    }

    ///@}

}; // class ReynoldsStressTensor
} // namespace Kratos

#endif // KRATOS_REYNOLDS_STRESS_TENSOR_H_INCLUDED defined