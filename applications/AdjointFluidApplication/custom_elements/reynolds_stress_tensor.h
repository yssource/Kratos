//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
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

    constexpr static double INV_TDIM = 1.0 / 3.0;

    ReynoldsStressTensor(Vector& rCoefficients,
                         Geometry<Node<3>>& rGeometry,
                         double TurbulentKineticEnergy,
                         double TurbulentKinematicViscosity,
                         double Density)
        : mCoefficients(rCoefficients), mGeometry(rGeometry)
    {
        mTurbulentKineticEnergy = TurbulentKineticEnergy;
        mTurbulentKinematicViscosity = TurbulentKinematicViscosity;
        mDensity = Density;

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            array_1d<double, 3>& rVel = mGeometry[i].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int j = 0; j < TDim; ++j)
                mVelocity(i, j) = rVel[j];
        }

        array_1d<double, TNumNodes> N;
        GeometryUtils::CalculateGeometryData(mGeometry, mIntegrationMatrix, N, mGaussPointWeight);

        KRATOS_ERROR_IF(mCoefficients.size() != 9)
            << "Non-linear eddy viscosity Reynold's stress model needs to have "
               "9 coefficients in the variable "
               "REYNOLDS_STRESS_MODEL_COEFFICIENTS. [Currently it has "
            << mCoefficients.size() << " coefficients.]";
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
                CalculateSymmetricVelocityGradientPrimalDerivative(
                    Deriv, SymmetricVelGradDerivative);

                const double SymmetricVelGradDerivTrace =
                    CalculateMatrixTrace(SymmetricVelGradDerivative);

                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int i = 0; i < TDim; ++i)
                    {
                        double valij = 0;
                        for (unsigned int l = 0; l < TDim; ++l)
                            valij += mTurbulentKinematicViscosity *
                                     SymmetricVelGradDerivative(l, i) *
                                     mIntegrationMatrix(a, l);

                        valij -= mTurbulentKinematicViscosity * INV_TDIM *
                                 mIntegrationMatrix(a, i) * SymmetricVelGradDerivTrace;

                        rLeftHandSideMatrix(a * TBlockSize + i, b * TBlockSize + j) +=
                            mDensity * valij * mGaussPointWeight;
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
                    rRightHandSideVector(a * TBlockSize + i) -=
                        mDensity * mIntegrationMatrix(a, l) *
                        reynolds_stress_tensor(l, i) * mGaussPointWeight;
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
                Deriv.NodeIndex = c;
                Deriv.Direction = k;
                this->CalculateReynoldsStressTensorPrimalDerivative(
                    Deriv, ReynoldsStressTensorDerivative);
                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int i = 0; i < TDim; ++i)
                    {
                        double valik = 0.0;
                        for (unsigned int l = 0; l < TDim; ++l)
                        {
                            valik += mIntegrationMatrix(a, l) *
                                     ReynoldsStressTensorDerivative(l, i);
                        }

                        rLeftHandSideMatrix(a * TBlockSize + i, c * TBlockSize + k) +=
                            mDensity * valik * mGaussPointWeight;
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
        Matrix DN_DX_Deriv(TNumNodes, TDim);
        double DetJ0_deriv;

        this->CalculateReynoldsStressTensor(ReynoldsStressTensorValue);

        for (unsigned int r = 0; r < TNumNodes; ++r)
        {
            for (unsigned int m = 0; m < TDim; ++m)
            {
                Deriv.NodeIndex = r;
                Deriv.Direction = m;
                this->CalculateReynoldsStressTensorShapeDerivative(
                    Deriv, ReynoldsStressTensorDerivative);
                this->CalculateShapeFunctionDerivativeShapeDerivatives(
                    0, Deriv, DetJ0_deriv, DN_DX_Deriv);

                for (unsigned int a = 0; a < TNumNodes; ++a)
                {
                    for (unsigned int i = 0; i < TDim; ++i)
                    {
                        double valim = 0.0;
                        for (unsigned int l = 0; l < TDim; ++l)
                        {
                            valim += mIntegrationMatrix(a, l) *
                                     ReynoldsStressTensorValue(l, i) * DetJ0_deriv * 0.5;
                            valim += mIntegrationMatrix(a, l) *
                                     ReynoldsStressTensorDerivative(l, i) * mGaussPointWeight;
                            valim += DN_DX_Deriv(a, l) *
                                     ReynoldsStressTensorValue(l, i) * mGaussPointWeight;
                        }

                        rShapeDerivativesMatrix(r * TDim + m,
                                                a * TBlockSize + i) -= mDensity * valim;
                    }
                }
            }
        }
    }

    void CalculateReynoldsStressTensor(BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        AddDeviatoricReynoldsStressTensorLinearPart(rOutput);
        // AddDeviatoricReynoldsStressTensorNonLinearPart(rOutput);
        // AddHydrostaticReynoldsStressTensor(rOutput);
    }

    void CalculateReynoldsStressTensorPrimalDerivative(const ShapeParameter& Deriv,
                                                       BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        AddDeviatoricReynoldsStressTensorLinearPartPrimalDerivative(Deriv, rOutput);
        AddDeviatoricReynoldsStressTensorNonLinearPartPrimalDerivative(
            Deriv, rOutput);
        AddHydrostaticReynoldsStressTensorPrimalDerivative(Deriv, rOutput);
    }

    void CalculateReynoldsStressTensorShapeDerivative(const ShapeParameter& Deriv,
                                                      BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        rOutput.clear();
        AddDeviatoricReynoldsStressTensorLinearPartShapeDerivative(Deriv, rOutput);
        AddDeviatoricReynoldsStressTensorNonLinearPartShapeDerivative(
            Deriv, rOutput);
        AddHydrostaticReynoldsStressTensorShapeDerivative(Deriv, rOutput);
    }

private:
    ///@name Member Variables
    ///@{
    BoundedMatrix<double, TNumNodes, TDim> mVelocity;
    BoundedMatrix<double, TNumNodes, TDim> mIntegrationMatrix;
    Vector& mCoefficients;
    Geometry<Node<3>>& mGeometry;
    double mDensity;
    double mGaussPointWeight;
    double mTurbulentKineticEnergy;
    double mTurbulentKinematicViscosity;

    double CalculateMatrixTrace(const Matrix& rMatrix)
    {
        double result = 0.0;
        for (unsigned int i = 0; i < rMatrix.size1(); i++)
            result += rMatrix(i, i);
        return result;
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

        const auto& rGeom = mGeometry;

        Geometry<Point>::JacobiansType J;
        rGeom.Jacobian(J, GeometryData::GI_GAUSS_1);

        Geometry<Point>::ShapeFunctionsGradientsType DN_De;
        DN_De = rGeom.ShapeFunctionsLocalGradients(GeometryData::GI_GAUSS_1);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

        const Matrix& rJ = J[0];
        const Matrix& rDN_De = DN_De[0];

        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        geom_sensitivity.CalculateSensitivity(Deriv, rDetJ0_deriv, rDN_DX0_deriv);

        KRATOS_CATCH("")
    }

    void CalculateSymmetricVelocityGradientShapeDerivative(
        const ShapeParameter Deriv, BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        rOutput.clear();

        Matrix DN_DX0_deriv(TNumNodes, TDim);
        double DetJ0_deriv;
        CalculateShapeFunctionDerivativeShapeDerivatives(0, Deriv, DetJ0_deriv, DN_DX0_deriv);

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
            {
                for (unsigned int d = 0; d < TNumNodes; d++)
                {
                    rOutput(i, j) += mVelocity(d, i) * DN_DX0_deriv(d, j);
                    rOutput(i, j) += mVelocity(d, j) * DN_DX0_deriv(d, i);
                }
            }

        rOutput *= 0.5;

        KRATOS_CATCH("")
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

    void AddDeviatoricReynoldsStressTensorNonLinearPart(
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        BoundedMatrix<double, 3, 3> reynolds_stress_tensor;

        for (unsigned int i = 0; i < 3; i++)
            for (unsigned int j = 0; j < 3; j++)
                reynolds_stress_tensor(i, j) = mCoefficients[i * 3 + j];

        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                rOutput(i, j) += reynolds_stress_tensor(i, j);

        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorLinearPart(
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);

        BoundedMatrix<double, TDim, TDim> SymmetricVelGrad;
        CalculateSymmetricVelocityGradient(this->mVelocity, SymmetricVelGrad);

        rOutput += mTurbulentKinematicViscosity *
                   (SymmetricVelGrad -
                    INV_TDIM * CalculateMatrixTrace(SymmetricVelGrad) * mIdentity);

        KRATOS_CATCH("")
    }

    void AddHydrostaticReynoldsStressTensor(BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);
        rOutput -= 2 * INV_TDIM * mTurbulentKineticEnergy * mIdentity;

        KRATOS_CATCH("")
    }

    void AddHydrostaticReynoldsStressTensorPrimalDerivative(
        const ShapeParameter& Deriv,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY
        // Nothing to add, since hydro static part of the reynolds stress tensor is constant
        KRATOS_CATCH("")
    }

    void AddHydrostaticReynoldsStressTensorShapeDerivative(
        const ShapeParameter& Deriv,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY
        // Nothing to add, since hydro static part of the reynolds stress tensor is constant
        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorNonLinearPartPrimalDerivative(
        const ShapeParameter& Deriv,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY
        // Nothing to add, since non-linear part of the reynolds stress tensor is constant
        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorNonLinearPartShapeDerivative(
        const ShapeParameter& Deriv,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY
        // Nothing to add, since non-linear part of the reynolds stress tensor is constant
        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorLinearPartPrimalDerivative(
        const ShapeParameter& Deriv,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);

        BoundedMatrix<double, TDim, TDim> SymmetricVelGradPrimalDerivative;
        CalculateSymmetricVelocityGradientPrimalDerivative(
            Deriv, SymmetricVelGradPrimalDerivative);

        rOutput += mTurbulentKinematicViscosity *
                   (SymmetricVelGradPrimalDerivative -
                    INV_TDIM * CalculateMatrixTrace(SymmetricVelGradPrimalDerivative) * mIdentity);

        KRATOS_CATCH("")
    }

    void AddDeviatoricReynoldsStressTensorLinearPartShapeDerivative(
        const ShapeParameter& Deriv,
        BoundedMatrix<double, TDim, TDim>& rOutput)
    {
        KRATOS_TRY

        IdentityMatrix mIdentity(TDim);

        BoundedMatrix<double, TDim, TDim> SymmetricVelGradShapeDerivative;
        CalculateSymmetricVelocityGradientShapeDerivative(Deriv, SymmetricVelGradShapeDerivative);

        rOutput += mTurbulentKinematicViscosity *
                   (SymmetricVelGradShapeDerivative -
                    INV_TDIM * CalculateMatrixTrace(SymmetricVelGradShapeDerivative) * mIdentity);

        KRATOS_CATCH("")
    }

    ///@}

}; // class ReynoldsStressTensor
} // namespace Kratos

#endif // KRATOS_REYNOLDS_STRESS_TENSOR_H_INCLUDED defined