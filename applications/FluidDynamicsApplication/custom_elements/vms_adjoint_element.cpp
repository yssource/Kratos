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
#include "utilities/adjoint_extensions.h"
#include "utilities/geometry_utilities.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "vms_adjoint_element.h"

namespace Kratos
{
///@name generic implementations

template <unsigned int TDim>
VMSAdjointElement<TDim>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension() + 1);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Y, Step);
    if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
    {
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Z, Step);
    }
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension() + 1);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Y, Step);
    if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
    {
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Z, Step);
    }
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension() + 1);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Y, Step);
    if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
    {
        rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Z, Step);
    }
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
}

template <unsigned int TDim>
VMSAdjointElement<TDim>::VMSAdjointElement(IndexType NewId) : Element(NewId)
{
}

template <unsigned int TDim>
VMSAdjointElement<TDim>::VMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

template <unsigned int TDim>
VMSAdjointElement<TDim>::VMSAdjointElement(IndexType NewId,
                                           GeometryType::Pointer pGeometry,
                                           PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim>
VMSAdjointElement<TDim>::~VMSAdjointElement()
{
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::Initialize()
{
    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
}

template <unsigned int TDim>
Element::Pointer VMSAdjointElement<TDim>::Create(IndexType NewId,
                                                 NodesArrayType const& ThisNodes,
                                                 PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Element::Pointer(new VMSAdjointElement<TDim>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));

    KRATOS_CATCH("")
}

template <unsigned int TDim>
Element::Pointer VMSAdjointElement<TDim>::Create(IndexType NewId,
                                                 GeometryType::Pointer pGeom,
                                                 PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Element::Pointer(new VMSAdjointElement<TDim>(NewId, pGeom, pProperties));

    KRATOS_CATCH("")
}

template <unsigned int TDim>
int VMSAdjointElement<TDim>::Check(const ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    // Check the element id and geometry.
    ProcessInfo UnusedProcessInfo;
    int ReturnValue = Element::Check(UnusedProcessInfo);

    // Check if adjoint and fluid variables are defined.
    if (ADJOINT_FLUID_VECTOR_1.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ADJOINT_FLUID_VECTOR_1 Key is 0. "
                           "Check if the application was correctly registered.",
                           "");
    if (ADJOINT_FLUID_VECTOR_2.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ADJOINT_FLUID_VECTOR_2 Key is 0. "
                           "Check if the application was correctly registered.",
                           "");
    if (ADJOINT_FLUID_VECTOR_3.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ADJOINT_FLUID_VECTOR_3 Key is 0. "
                           "Check if the application was correctly registered.",
                           "");
    if (ADJOINT_FLUID_SCALAR_1.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ADJOINT_FLUID_SCALAR_1 Key is 0. "
                           "Check if the application was correctly registered.",
                           "");
    if (VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "VELOCITY Key is 0. "
                           "Check if the application was correctly registered.",
                           "");
    if (ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ACCELERATION Key is 0. "
                           "Check if the application was correctly registered.",
                           "");
    if (PRESSURE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "PRESSURE Key is 0. "
                           "Check if the application was correctly registered.",
                           "");

    // Check if the nodes have adjoint and fluid variables and adjoint dofs.
    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_1) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing ADJOINT_FLUID_VECTOR_1 variable on "
                               "solution step data for node ",
                               this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_2) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing ADJOINT_FLUID_VECTOR_2 variable on "
                               "solution step data for node ",
                               this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_3) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing ADJOINT_FLUID_VECTOR_3 variable on "
                               "solution step data for node ",
                               this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_SCALAR_1) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing ADJOINT_FLUID_SCALAR_1 variable on "
                               "solution step data for node ",
                               this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(
                std::invalid_argument,
                "missing VELOCITY variable on solution step data for node ",
                this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(ACCELERATION) == false)
            KRATOS_THROW_ERROR(
                std::invalid_argument,
                "missing ACCELERATION variable on solution step data for node ",
                this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_THROW_ERROR(
                std::invalid_argument,
                "missing PRESSURE variable on solution step data for node ",
                this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_VECTOR_1_X) == false ||
            this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_VECTOR_1_Y) == false ||
            this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_VECTOR_1_Z) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing ADJOINT_FLUID_VECTOR_1 component "
                               "degree of freedom on node ",
                               this->GetGeometry()[iNode].Id());
        if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_SCALAR_1) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing ADJOINT_FLUID_SCALAR_1 component "
                               "degree of freedom on node ",
                               this->GetGeometry()[iNode].Id());
    }

    return ReturnValue;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::GetValuesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TFluidLocalSize)
        rValues.resize(TFluidLocalSize, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        const array_1d<double, 3>& rVel =
            rGeom[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, Step);
        for (IndexType d = 0; d < TDim; d++)
            rValues[LocalIndex++] = rVel[d];
        rValues[LocalIndex++] =
            rGeom[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
    }
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::GetFirstDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TFluidLocalSize)
        rValues.resize(TFluidLocalSize, false);

    rValues.clear();
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::GetSecondDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TFluidLocalSize)
        rValues.resize(TFluidLocalSize, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        const array_1d<double, 3>& rAccel =
            rGeom[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, Step);
        for (IndexType d = 0; d < TDim; d++)
            rValues[LocalIndex++] = rAccel[d];
        rValues[LocalIndex++] = 0.0; // pressure dof
    }
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                   VectorType& rRightHandSideVector,
                                                   ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR(std::runtime_error, "this function is not implemented.",
                       "")

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                    ProcessInfo& /*rCurrentProcessInfo*/)
{
    if (rLeftHandSideMatrix.size1() != TFluidLocalSize ||
        rLeftHandSideMatrix.size2() != TFluidLocalSize)
        rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                     ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR(std::runtime_error, "this function is not implemented.",
                       "")

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                           ProcessInfo& rCurrentProcessInfo)
{
    this->CalculatePrimalGradientOfVMSSteadyTerm(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->AddPrimalGradientOfVMSMassTerm(rLeftHandSideMatrix, ACCELERATION,
                                         -1.0, rCurrentProcessInfo);
    rLeftHandSideMatrix = trans(rLeftHandSideMatrix); // transpose
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                            ProcessInfo& rCurrentProcessInfo)
{
    this->CalculateVMSMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    rLeftHandSideMatrix = -trans(rLeftHandSideMatrix); // transpose
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                  ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR(std::runtime_error, "this function is not implemented.",
                       "")

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                                     ProcessInfo& /*rCurrentProcessInfo*/)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR(std::runtime_error, "this function is not implemented.",
                       "")

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY)
    {
        this->CalculateShapeGradientOfVMSSteadyTerm(rOutput, rCurrentProcessInfo);
        this->AddShapeGradientOfVMSMassTerm(rOutput, ACCELERATION, -1.0, rCurrentProcessInfo);
    }
    else
    {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
std::string VMSAdjointElement<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "VMSAdjointElement" << this->GetGeometry().WorkingSpaceDimension()
           << "D #" << this->Id();
    return buffer.str();
}
template <unsigned int TDim>
void VMSAdjointElement<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "VMSAdjointElement" << this->GetGeometry().WorkingSpaceDimension()
             << "D #" << this->Id() << std::endl;
    rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber() << std::endl;
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::PrintData(std::ostream& rOStream) const
{
    this->PrintInfo(rOStream);
    rOStream << "Geometry Data: " << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateVMSMassMatrix(MatrixType& rMassMatrix,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rMassMatrix.size1() != TFluidLocalSize || rMassMatrix.size2() != TFluidLocalSize)
        rMassMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);

    rMassMatrix.clear();

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density, DENSITY, N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity, VISCOSITY, N);
    Viscosity *= Density;

    // u
    array_1d<double, TDim> Velocity;
    this->EvaluateInPoint(Velocity, VELOCITY, N);

    // u * Grad(N)
    array_1d<double, TNumNodes> DensityVelGradN;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        DensityVelGradN[i] = 0.0;
        for (IndexType d = 0; d < TDim; d++)
            DensityVelGradN[i] += Density * DN_DX(i, d) * Velocity[d];
    }

    // Stabilization parameters
    double VelNorm = 0.0;
    for (IndexType d = 0; d < TDim; d++)
        VelNorm += Velocity[d] * Velocity[d];
    VelNorm = std::sqrt(VelNorm);
    const double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize,
                                           Density, Viscosity, rCurrentProcessInfo);

    // Lumped mass
    const double LumpedMass = Density * Volume / static_cast<double>(TNumNodes);
    IndexType DofIndex = 0;
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType d = 0; d < TDim; ++d)
        {
            rMassMatrix(DofIndex, DofIndex) += LumpedMass;
            ++DofIndex;
        }
        ++DofIndex; // Skip pressure Dof
    }

    // Stabilization, convection-acceleration
    IndexType FirstRow(0), FirstCol(0);
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const double diag = DensityVelGradN[i] * TauOne * Density * N[j];

            for (IndexType d = 0; d < TDim; ++d)
            {
                rMassMatrix(FirstRow + d, FirstCol + d) += Volume * diag;
                rMassMatrix(FirstRow + TDim, FirstCol + d) +=
                    Volume * DN_DX(i, d) * TauOne * Density * N[j];
            }

            FirstCol += TBlockSize;
        } // Node block columns

        FirstRow += TBlockSize;
        FirstCol = 0;
    } // Node block rows

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::AddPrimalGradientOfVMSMassTerm(
    MatrixType& rOutputMatrix,
    const Variable<array_1d<double, 3>>& rVariable,
    double alpha,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rOutputMatrix.size1() != TFluidLocalSize)
    {
        KRATOS_THROW_ERROR(
            std::runtime_error,
            "invalid matrix size detected. rOutputMatrix.size1() = ",
            rOutputMatrix.size1());
    }

    if (rOutputMatrix.size2() != TFluidLocalSize)
    {
        KRATOS_THROW_ERROR(
            std::runtime_error,
            "invalid matrix size detected. rOutputMatrix.size2() = ",
            rOutputMatrix.size2());
    }

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density, DENSITY, N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity, VISCOSITY, N);
    Viscosity *= Density;

    // u
    array_1d<double, TDim> Velocity;
    this->EvaluateInPoint(Velocity, VELOCITY, N);

    // u * Grad(N)
    array_1d<double, TNumNodes> DensityVelGradN;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        DensityVelGradN[i] = 0.0;
        for (IndexType d = 0; d < TDim; d++)
            DensityVelGradN[i] += Density * DN_DX(i, d) * Velocity[d];
    }

    // Stabilization parameters
    double VelNorm = 0.0;
    for (IndexType d = 0; d < TDim; d++)
        VelNorm += Velocity[d] * Velocity[d];
    VelNorm = std::sqrt(VelNorm);
    const double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize,
                                           Density, Viscosity, rCurrentProcessInfo);

    // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
    // depend on the definitions of TauOne and TauTwo and should be consistent
    // with the fluid element used to solve for VELOCITY and PRESSURE.
    BoundedMatrix<double, TNumNodes, TDim> TauOneDeriv;
    BoundedMatrix<double, TNumNodes, TDim> TauTwoDeriv;

    if (VelNorm > 0.0)
    {
        const double CoefOne = -2.0 * Density * TauOne * TauOne / (ElemSize * VelNorm);
        const double CoefTwo = 0.5 * Density * ElemSize / VelNorm;

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                TauOneDeriv(i, d) = CoefOne * N[i] * Velocity[d];
                TauTwoDeriv(i, d) = CoefTwo * N[i] * Velocity[d];
            }
        }
    }

    // rVariable (x)
    array_1d<double, TDim> X;
    this->EvaluateInPoint(X, rVariable, N);

    // x * Grad(N)
    array_1d<double, TNumNodes> DensityXGradN;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        DensityXGradN[i] = 0.0;
        for (IndexType d = 0; d < TDim; d++)
            DensityXGradN[i] += Density * DN_DX(i, d) * X[d];
    }

    // Primal gradient of (lumped) VMS mass matrix multiplied with vector
    IndexType FirstRow(0), FirstCol(0);
    // Loop over nodes
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            for (IndexType m = 0; m < TDim; ++m)
            {
                for (IndexType n = 0; n < TDim; ++n)
                {
                    double valmn = 0.0;

                    valmn += DensityVelGradN[i] * TauOneDeriv(j, n) * Density * X[m];

                    valmn += Density * N[j] * DN_DX(i, n) * TauOne * Density * X[m];

                    rOutputMatrix(FirstRow + m, FirstCol + n) += alpha * Volume * valmn;
                }

                rOutputMatrix(FirstRow + TDim, FirstCol + m) +=
                    alpha * Volume * DensityXGradN[i] * TauOneDeriv(j, m);
            }

            FirstCol += TBlockSize;
        } // Node block columns

        FirstRow += TBlockSize;
        FirstCol = 0;
    } // Node block rows

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::AddShapeGradientOfVMSMassTerm(
    MatrixType& rOutputMatrix,
    const Variable<array_1d<double, 3>>& rVariable,
    double alpha,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rOutputMatrix.size1() != TCoordLocalSize)
    {
        KRATOS_THROW_ERROR(
            std::runtime_error,
            "invalid matrix size detected. rOutputMatrix.size1() = ",
            rOutputMatrix.size1());
    }

    if (rOutputMatrix.size2() != TFluidLocalSize)
    {
        KRATOS_THROW_ERROR(
            std::runtime_error,
            "invalid matrix size detected. rOutputMatrix.size2() = ",
            rOutputMatrix.size2());
    }

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density, DENSITY, N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity, VISCOSITY, N);
    Viscosity *= Density;

    // u
    array_1d<double, TDim> Velocity;
    this->EvaluateInPoint(Velocity, VELOCITY, N);

    // u * Grad(N)
    array_1d<double, TNumNodes> DensityVelGradN;
    noalias(DensityVelGradN) = Density * prod(DN_DX, Velocity);

    // Det(J)
    const double InvDetJ = 1.0 / this->GetGeometry().DeterminantOfJacobian(0);
    array_1d<double, TCoordLocalSize> DetJDerivatives;
    this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

    // Stabilization parameters TauOne, TauTwo
    double VelNorm = norm_2(Velocity);
    double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize,
                                           Density, Viscosity, rCurrentProcessInfo);

    // Vector values
    array_1d<double, TFluidLocalSize> X;
    IndexType DofIndex = 0;
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        array_1d<double, 3>& rValue =
            this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);
        for (IndexType d = 0; d < TDim; ++d)
        {
            X[DofIndex++] = rValue[d];
        }
        X[DofIndex++] = 0.0; // pressure dof
    }

    array_1d<double, TFluidLocalSize> Derivative;

    // We compute the derivative w.r.t each coordinate of each node and
    // assign it to the corresponding row of the shape derivatives matrix.
    for (IndexType iCoord = 0; iCoord < TCoordLocalSize; ++iCoord)
    {
        // Det(J)'
        double DetJDeriv = DetJDerivatives[iCoord];

        // DN_DX'
        BoundedMatrix<double, TNumNodes, TDim> DN_DX_Deriv;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                DN_DX_Deriv(i, d) = -DN_DX(iCoord / TDim, d) * DN_DX(i, iCoord % TDim);
            }
        }

        // Volume'
        double VolumeDeriv = Volume * InvDetJ * DetJDeriv;

        // u * Grad(N)'
        array_1d<double, TNumNodes> DensityVelGradNDeriv;
        noalias(DensityVelGradNDeriv) = Density * prod(DN_DX_Deriv, Velocity);

        // TauOne', TauTwo'
        double TauOneDeriv, TauTwoDeriv;
        this->CalculateStabilizationParametersDerivative(
            TauOneDeriv, TauTwoDeriv, TauOne, TauTwo, VelNorm, ElemSize,
            Density, Viscosity, DetJDeriv);

        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> LHS;
        array_1d<double, TFluidLocalSize> RHS;
        for (IndexType i = 0; i < TFluidLocalSize; i++)
        {
            RHS[i] = 0.0;
            for (IndexType j = 0; j < TFluidLocalSize; j++)
                LHS(i, j) = 0.0;
        }

        // The usual lumped mass matrix
        double LumpedMassDeriv = Density * VolumeDeriv / static_cast<double>(TNumNodes);
        IndexType DofIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                LHS(DofIndex, DofIndex) += LumpedMassDeriv;
                ++DofIndex;
            }
            ++DofIndex; // Skip pressure Dof
        }

        // Stabilization, convection-acceleration
        IndexType FirstRow(0), FirstCol(0);
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                double diag = 0.0;
                double ddiag = 0.0;

                diag += DensityVelGradN[i] * TauOne * Density * N[j];
                ddiag += DensityVelGradNDeriv[i] * TauOne * Density * N[j];
                ddiag += DensityVelGradN[i] * TauOneDeriv * Density * N[j];

                for (IndexType n = 0; n < TDim; ++n)
                {
                    double valn = DN_DX(i, n) * TauOne * Density * N[j];
                    double dvaln = 0.0;
                    dvaln += DN_DX_Deriv(i, n) * TauOne * Density * N[j];
                    dvaln += DN_DX(i, n) * TauOneDeriv * Density * N[j];

                    LHS(FirstRow + n, FirstCol + n) += VolumeDeriv * diag + Volume * ddiag;

                    LHS(FirstRow + TDim, FirstCol + n) +=
                        VolumeDeriv * valn + Volume * dvaln;
                }

                FirstCol += TBlockSize;
            } // Node block columns

            FirstRow += TBlockSize;
            FirstCol = 0;
        } // Node block rows

        // Assign the derivative w.r.t this coordinate to the
        // shape derivative mass matrix.
        noalias(Derivative) = prod(LHS, X);
        for (IndexType k = 0; k < TFluidLocalSize; ++k)
        {
            rOutputMatrix(iCoord, k) += alpha * Derivative[k];
        }
    }
    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculatePrimalGradientOfVMSSteadyTerm(MatrixType& rAdjointMatrix,
                                                                     const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rAdjointMatrix.size1() != TFluidLocalSize || rAdjointMatrix.size2() != TFluidLocalSize)
        rAdjointMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);

    rAdjointMatrix.clear();

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density, DENSITY, N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity, VISCOSITY, N);
    Viscosity *= Density;

    // u
    array_1d<double, TDim> Velocity;
    this->EvaluateInPoint(Velocity, VELOCITY, N);

    // u * Grad(N)
    array_1d<double, TNumNodes> DensityVelGradN;
    noalias(DensityVelGradN) = Density * prod(DN_DX, Velocity);

    // Grad(u)
    BoundedMatrix<double, TDim, TDim> DensityGradVel;
    this->CalculateVelocityGradient(DensityGradVel, DN_DX);

    // Div(u)
    double DivVel = 0.0;
    for (IndexType d = 0; d < TDim; ++d)
        DivVel += DensityGradVel(d, d);

    DensityGradVel *= Density;

    // Grad(p)
    array_1d<double, TDim> GradP;
    this->CalculatePressureGradient(GradP, DN_DX);

    // ( Grad(u) * Grad(N) )^T
    BoundedMatrix<double, TNumNodes, TDim> DN_DX_DensityGradVel;
    noalias(DN_DX_DensityGradVel) = prod(DN_DX, DensityGradVel);

    // ( u * Grad(u) * Grad(N) )^T
    array_1d<double, TNumNodes> DN_DX_DensityGradVel_Vel;
    noalias(DN_DX_DensityGradVel_Vel) = prod(DN_DX_DensityGradVel, Velocity);

    // u * Grad(u)
    array_1d<double, TDim> DensityGradVel_Vel;
    noalias(DensityGradVel_Vel) = prod(DensityGradVel, Velocity);

    // Grad(N)^T * Grad(p)
    array_1d<double, TNumNodes> DN_DX_GradP;
    noalias(DN_DX_GradP) = prod(DN_DX, GradP);

    // Grad(N)^T * BodyForce
    array_1d<double, TDim> BodyForce;
    array_1d<double, TNumNodes> DN_DX_BodyForce;
    this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
    BodyForce *= Density;
    noalias(DN_DX_BodyForce) = prod(DN_DX, BodyForce);

    // Stabilization parameters TauOne, TauTwo
    double VelNorm = norm_2(Velocity);
    double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize,
                                           Density, Viscosity, rCurrentProcessInfo);

    // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
    // depend on the definitions of TauOne and TauTwo and should be
    // consistent with the fluid element used to solve for VELOCITY and
    // PRESSURE.
    BoundedMatrix<double, TNumNodes, TDim> TauOneDeriv;
    BoundedMatrix<double, TNumNodes, TDim> TauTwoDeriv;

    if (VelNorm > 0.0)
    {
        double CoefOne = -2.0 * Density * TauOne * TauOne / (ElemSize * VelNorm);
        double CoefTwo = 0.5 * Density * ElemSize / VelNorm;

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                TauOneDeriv(i, d) = CoefOne * N[i] * Velocity[d];
                TauTwoDeriv(i, d) = CoefTwo * N[i] * Velocity[d];
            }
        }
    }

    // Here, -(\partial R / \partial W) is calculated. This is the discrete
    // derivative of the fluid residual w.r.t the fluid variables and therefore
    // includes many of the terms defined in the fluid element. Neglecting the
    // transient terms of the fluid element, this matrix is identical to the
    // Jacobian of the fluid residual used for Newton-Raphson iterations. The
    // matrix is transposed at the end to get the adjoint system matrix.

    IndexType FirstRow(0), FirstCol(0);
    // Loop over nodes
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            double diag = 0.0;

            // Convective term, v * (u * Grad(u))
            diag += N[i] * DensityVelGradN[j];

            // Stabilization, lsq convection
            // (u * Grad(v)) * TauOne * (u * Grad(u))
            diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];

            for (IndexType m = 0; m < TDim; ++m)
            {
                for (IndexType n = 0; n < TDim; ++n)
                {
                    double valmn = 0.0;

                    // Convective term, v * (u * Grad(u))
                    valmn += N[i] * N[j] * DensityGradVel(m, n);

                    // Stabilization, lsq convection
                    // (u * Grad(v)) * TauOne * (u * Grad(u))
                    valmn += DensityVelGradN[i] * TauOne * N[j] * DensityGradVel(m, n);
                    valmn += DensityVelGradN[i] * TauOneDeriv(j, n) *
                             DensityGradVel_Vel[m];
                    valmn += Density * N[j] * DN_DX(i, n) * TauOne * DensityGradVel_Vel[m];

                    // Stabilization, lsq divergence
                    // Div(v) * TauTwo * Div(u)
                    valmn += DN_DX(i, m) * TauTwo * DN_DX(j, n);
                    valmn += DN_DX(i, m) * TauTwoDeriv(j, n) * DivVel;

                    // Stabilization, convection-pressure
                    // (u * Grad(v)) * TauOne * Grad(p)
                    valmn += TauOneDeriv(j, n) * DensityVelGradN[i] * GradP[m];
                    valmn += Density * TauOne * N[j] * DN_DX(i, n) * GradP[m];

                    // Stabilization, convection-BodyForce
                    // (u * Grad(v)) * TauOne * f
                    valmn -= N[j] * DN_DX(i, n) * TauOne * Density * BodyForce[m];
                    valmn -= DensityVelGradN[i] * TauOneDeriv(j, n) * BodyForce[m];

                    rAdjointMatrix(FirstRow + m, FirstCol + n) += Volume * valmn;
                }

                rAdjointMatrix(FirstRow + m, FirstCol + m) += Volume * diag;

                double valmp = 0.0;
                double valpn = 0.0;

                // Pressure term
                // Div(v) * p
                valmp -= DN_DX(i, m) * N[j];

                // Stabilization, convection-pressure
                // (u * Grad(v)) * TauOne * Grad(p)
                valmp += TauOne * DensityVelGradN[i] * DN_DX(j, m);

                // Divergence term
                // q * Div(u)
                valpn += N[i] * DN_DX(j, m);

                // Stabilization, lsq pressure
                // TauOne * Grad(q) * Grad(p)
                valpn += DN_DX_GradP[i] * TauOneDeriv(j, m);

                // Stabilization, pressure-convection
                // Grad(q) * TauOne * (u * Grad(u))
                valpn += DN_DX(i, m) * TauOne * DensityVelGradN[j];
                valpn += DN_DX_DensityGradVel(i, m) * TauOne * N[j];
                valpn += DN_DX_DensityGradVel_Vel[i] * TauOneDeriv(j, m);

                // Stabilization, pressure-BodyForce
                // Grad(q) * TauOne * f
                valpn -= DN_DX_BodyForce[i] * TauOneDeriv(j, m);

                rAdjointMatrix(FirstRow + m, FirstCol + TDim) += Volume * valmp;
                rAdjointMatrix(FirstRow + TDim, FirstCol + m) += Volume * valpn;
            }

            // Stabilization, lsq pressure
            // TauOne * Grad(q) * Grad(p)
            double valpp = 0.0;
            for (IndexType d = 0; d < TDim; ++d)
            {
                valpp += DN_DX(i, d) * DN_DX(j, d);
            }
            valpp *= TauOne;

            rAdjointMatrix(FirstRow + TDim, FirstCol + TDim) += Volume * valpp;

            FirstCol += TBlockSize;
        } // Node block columns

        FirstRow += TBlockSize;
        FirstCol = 0;
    } // Node block rows

    // Viscous term
    this->AddViscousTerm(rAdjointMatrix, DN_DX, Viscosity * Volume);

    // change the sign for consistency with definition
    noalias(rAdjointMatrix) = -rAdjointMatrix;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateShapeGradientOfVMSSteadyTerm(
    MatrixType& rShapeDerivativesMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rShapeDerivativesMatrix.size1() != TCoordLocalSize ||
        rShapeDerivativesMatrix.size2() != TFluidLocalSize)
    {
        rShapeDerivativesMatrix.resize(TCoordLocalSize, TFluidLocalSize, false);
    }

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType DN_DX;
    array_1d<double, TNumNodes> N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density, DENSITY, N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity, VISCOSITY, N);
    Viscosity *= Density;

    // u
    array_1d<double, TDim> Velocity;
    this->EvaluateInPoint(Velocity, VELOCITY, N);

    // u * Grad(N)
    array_1d<double, TNumNodes> DensityVelGradN;
    noalias(DensityVelGradN) = Density * prod(DN_DX, Velocity);

    // Det(J)
    const double InvDetJ = 1.0 / this->GetGeometry().DeterminantOfJacobian(0);
    array_1d<double, TCoordLocalSize> DetJDerivatives;
    this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

    // Stabilization parameters TauOne, TauTwo
    double VelNorm = norm_2(Velocity);
    double ElemSize = this->CalculateElementSize(Volume);
    double TauOne, TauTwo;
    this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize,
                                           Density, Viscosity, rCurrentProcessInfo);

    // External body force
    array_1d<double, TDim> BodyForce;
    this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
    BodyForce *= Density;

    array_1d<double, TFluidLocalSize> FluidValues;

    IndexType DofIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; iNode++)
    {
        array_1d<double, 3>& rVelocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
        for (IndexType d = 0; d < TDim; d++)
            FluidValues[DofIndex++] = rVelocity[d];
        FluidValues[DofIndex++] =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE);
    }

    // We compute the derivative of the residual w.r.t each coordinate of
    // each node and assign it to the corresponding row of the shape
    // derivatives matrix.
    for (IndexType iCoord = 0; iCoord < TCoordLocalSize; ++iCoord)
    {
        // Det(J)'
        double DetJDeriv = DetJDerivatives[iCoord];

        // DN_DX'
        BoundedMatrix<double, TNumNodes, TDim> DN_DX_Deriv;
        for (IndexType i = 0; i < TNumNodes; ++i)
            for (IndexType d = 0; d < TDim; ++d)
                DN_DX_Deriv(i, d) = -DN_DX(iCoord / TDim, d) * DN_DX(i, iCoord % TDim);

        // Volume'
        double VolumeDeriv = Volume * InvDetJ * DetJDeriv;

        // u * Grad(N)'
        array_1d<double, TNumNodes> DensityVelGradNDeriv;
        noalias(DensityVelGradNDeriv) = Density * prod(DN_DX_Deriv, Velocity);

        // TauOne', TauTwo'
        double TauOneDeriv, TauTwoDeriv;
        this->CalculateStabilizationParametersDerivative(
            TauOneDeriv, TauTwoDeriv, TauOne, TauTwo, VelNorm, ElemSize,
            Density, Viscosity, DetJDeriv);

        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> LHS;
        array_1d<double, TFluidLocalSize> RHS;
        for (IndexType i = 0; i < TFluidLocalSize; i++)
        {
            RHS[i] = 0.0;
            for (IndexType j = 0; j < TFluidLocalSize; j++)
                LHS(i, j) = 0.0;
        }

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                // Left-hand side matrix
                double diag = 0.0;
                double ddiag = 0.0;

                // Convective term, v * (u * Grad(u))
                diag += N[i] * DensityVelGradN[j];
                ddiag += N[i] * DensityVelGradNDeriv[j];

                // Stabilization, lsq convection
                // (u * Grad(v)) * TauOne * (u * Grad(u))
                diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];
                ddiag += DensityVelGradNDeriv[i] * TauOne * DensityVelGradN[j] +
                         DensityVelGradN[i] * TauOneDeriv * DensityVelGradN[j] +
                         DensityVelGradN[i] * TauOne * DensityVelGradNDeriv[j];

                for (IndexType m = 0; m < TDim; ++m)
                {
                    for (IndexType n = 0; n < TDim; ++n)
                    {
                        // Stabilization, lsq divergence
                        // Div(v) * TauTwo * Div(u)
                        double valmn = DN_DX(i, m) * TauTwo * DN_DX(j, n);
                        double dvalmn = DN_DX_Deriv(i, m) * TauTwo * DN_DX(j, n) +
                                        DN_DX(i, m) * TauTwoDeriv * DN_DX(j, n) +
                                        DN_DX(i, m) * TauTwo * DN_DX_Deriv(j, n);

                        LHS(i * TBlockSize + m, j * TBlockSize + n) +=
                            VolumeDeriv * valmn + Volume * dvalmn;
                    }
                    LHS(i * TBlockSize + m, j * TBlockSize + m) +=
                        VolumeDeriv * diag + Volume * ddiag;

                    double valmp = 0.0;
                    double dvalmp = 0.0;
                    // Pressure term
                    // Div(v) * p
                    valmp -= DN_DX(i, m) * N[j];
                    dvalmp -= DN_DX_Deriv(i, m) * N[j];

                    // Stabilization, convection-pressure
                    // (u * Grad(v)) * TauOne * Grad(p)
                    valmp += TauOne * DensityVelGradN[i] * DN_DX(j, m);
                    dvalmp += TauOneDeriv * DensityVelGradN[i] * DN_DX(j, m) +
                              TauOne * DensityVelGradNDeriv[i] * DN_DX(j, m) +
                              TauOne * DensityVelGradN[i] * DN_DX_Deriv(j, m);

                    double valpn = 0.0;
                    double dvalpn = 0.0;
                    // Divergence term
                    // q * Div(u)
                    valpn += N[i] * DN_DX(j, m);
                    dvalpn += N[i] * DN_DX_Deriv(j, m);

                    // Stabilization, pressure-convection
                    // Grad(q) * TauOne * (u * Grad(u))
                    valpn += TauOne * DensityVelGradN[j] * DN_DX(i, m);
                    dvalpn += TauOneDeriv * DensityVelGradN[j] * DN_DX(i, m) +
                              TauOne * DensityVelGradNDeriv[j] * DN_DX(i, m) +
                              TauOne * DensityVelGradN[j] * DN_DX_Deriv(i, m);

                    LHS(i * TBlockSize + m, j * TBlockSize + TDim) +=
                        VolumeDeriv * valmp + Volume * dvalmp;
                    LHS(i * TBlockSize + TDim, j * TBlockSize + m) +=
                        VolumeDeriv * valpn + Volume * dvalpn;
                }

                double valpp = 0.0;
                double dvalpp = 0.0;
                // Stabilization, lsq pressure
                // TauOne * Grad(q) * Grad(p)
                for (IndexType d = 0; d < TDim; ++d)
                {
                    valpp += DN_DX(i, d) * DN_DX(j, d) * TauOne;
                    dvalpp += DN_DX_Deriv(i, d) * DN_DX(j, d) * TauOne +
                              DN_DX(i, d) * DN_DX_Deriv(j, d) * TauOne +
                              DN_DX(i, d) * DN_DX(j, d) * TauOneDeriv;
                }

                LHS(i * TBlockSize + TDim, j * TBlockSize + TDim) +=
                    VolumeDeriv * valpp + Volume * dvalpp;
            } // Node block columns

            // Right-hand side vector
            double DN_DX_BodyForce = 0.0;
            double DN_DX_BodyForceDeriv = 0.0;
            for (IndexType d = 0; d < TDim; ++d)
            {
                DN_DX_BodyForce += DN_DX(i, d) * BodyForce[d];
                DN_DX_BodyForceDeriv += DN_DX_Deriv(i, d) * BodyForce[d];
            }

            for (IndexType m = 0; m < TDim; ++m)
            {
                double valm = 0.0;
                double dvalm = 0.0;

                // External body force
                valm += N[i] * BodyForce[m];

                // Stabilization, convection-BodyForce
                // (u * Grad(v)) * TauOne * f
                valm += TauOne * DensityVelGradN[i] * BodyForce[m];
                dvalm += TauOneDeriv * DensityVelGradN[i] * BodyForce[m] +
                         TauOne * DensityVelGradNDeriv[i] * BodyForce[m];

                RHS[i * TBlockSize + m] += VolumeDeriv * valm + Volume * dvalm;
            }

            double valp = TauOne * DN_DX_BodyForce;
            double dvalp = TauOneDeriv * DN_DX_BodyForce + TauOne * DN_DX_BodyForceDeriv;

            RHS[i * TBlockSize + TDim] += VolumeDeriv * valp + Volume * dvalp;
        } // Node block rows

        this->AddViscousTermDerivative(
            LHS, DN_DX, DN_DX_Deriv, Viscosity * Volume, Viscosity * VolumeDeriv);

        // Assign the derivative of the residual w.r.t this coordinate to
        // the shape derivatives matrix.
        array_1d<double, TFluidLocalSize> ResidualDerivative;
        noalias(ResidualDerivative) = RHS - prod(LHS, FluidValues);
        for (IndexType k = 0; k < TFluidLocalSize; ++k)
            rShapeDerivativesMatrix(iCoord, k) = ResidualDerivative[k];
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateVelocityGradient(BoundedMatrix<double, TDim, TDim>& rGradVel,
                                                        const ShapeFunctionDerivativesType& rDN_DX)
{
    GeometryType& rGeom = this->GetGeometry();
    // node 0
    const array_1d<double, 3>& rVel = rGeom[0].FastGetSolutionStepValue(VELOCITY, 0);
    for (IndexType m = 0; m < TDim; m++)
        for (IndexType n = 0; n < TDim; n++)
            rGradVel(m, n) = rDN_DX(0, n) * rVel[m];
    // node 1,2,...
    for (IndexType iNode = 1; iNode < TNumNodes; iNode++)
    {
        const array_1d<double, 3>& rVel =
            rGeom[iNode].FastGetSolutionStepValue(VELOCITY, 0);
        for (IndexType m = 0; m < TDim; m++)
            for (IndexType n = 0; n < TDim; n++)
                rGradVel(m, n) += rDN_DX(iNode, n) * rVel[m];
    }
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculatePressureGradient(array_1d<double, TDim>& rGradP,
                                                        const ShapeFunctionDerivativesType& rDN_DX)
{
    GeometryType& rGeom = this->GetGeometry();
    // node 0
    for (IndexType d = 0; d < TDim; d++)
        rGradP[d] = rDN_DX(0, d) * rGeom[0].FastGetSolutionStepValue(PRESSURE);
    // node 1,2,...
    for (IndexType iNode = 1; iNode < TNumNodes; iNode++)
        for (IndexType d = 0; d < TDim; ++d)
            rGradP[d] += rDN_DX(iNode, d) * rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::CalculateStabilizationParameters(double& rTauOne,
                                                               double& rTauTwo,
                                                               double VelNorm,
                                                               double ElemSize,
                                                               double Density,
                                                               double Viscosity,
                                                               const ProcessInfo& rCurrentProcessInfo)
{
    // assume DELTA_TIME < 0 !!!
    double tmp = -rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME];
    tmp += 2.0 * VelNorm / ElemSize;
    tmp *= Density;
    tmp += 4.0 * Viscosity / (ElemSize * ElemSize);
    rTauOne = 1.0 / tmp;
    rTauTwo = Viscosity + 0.5 * Density * ElemSize * VelNorm;
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::EvaluateInPoint(double& rResult,
                                              const Variable<double>& rVariable,
                                              const array_1d<double, TNumNodes>& rShapeFunc,
                                              IndexType step)
{
    GeometryType& rGeom = this->GetGeometry();
    rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(rVariable, step);
    for (IndexType iNode = 1; iNode < TNumNodes; ++iNode)
    {
        rResult += rShapeFunc[iNode] *
                   rGeom[iNode].FastGetSolutionStepValue(rVariable, step);
    }
}

template <unsigned int TDim>
void VMSAdjointElement<TDim>::EvaluateInPoint(array_1d<double, TDim>& rResult,
                                              const Variable<array_1d<double, 3>>& rVariable,
                                              const array_1d<double, TNumNodes>& rN,
                                              IndexType step)
{
    GeometryType& rGeom = this->GetGeometry();
    const array_1d<double, 3>& rNodalValue =
        rGeom[0].FastGetSolutionStepValue(rVariable, step);
    for (IndexType d = 0; d < TDim; d++)
        rResult[d] = rN[0] * rNodalValue[d];
    for (IndexType iNode = 1; iNode < TNumNodes; ++iNode)
    {
        const array_1d<double, 3>& rNodalValue =
            rGeom[iNode].FastGetSolutionStepValue(rVariable, step);
        for (IndexType d = 0; d < TDim; d++)
            rResult[d] += rN[iNode] * rNodalValue[d];
    }
}

///@}

///@name Specialized implementation for functions that depend on TDim
///@{

template <>
void VMSAdjointElement<2>::GetDofList(DofsVectorType& rElementalDofList,
                                      ProcessInfo& /*rCurrentProcessInfo*/)
{
    const unsigned int NumNodes(3), LocalSize(9);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_VECTOR_1_X);
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_VECTOR_1_Y);
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_SCALAR_1);
    }
}

template <>
void VMSAdjointElement<3>::GetDofList(DofsVectorType& rElementalDofList,
                                      ProcessInfo& /*rCurrentProcessInfo*/)
{
    const SizeType NumNodes(4), LocalSize(16);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    IndexType LocalIndex = 0;

    for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_VECTOR_1_X);
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_VECTOR_1_Y);
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_VECTOR_1_Z);
        rElementalDofList[LocalIndex++] =
            this->GetGeometry()[iNode].pGetDof(ADJOINT_FLUID_SCALAR_1);
    }
}

template <>
void VMSAdjointElement<2>::EquationIdVector(EquationIdVectorType& rResult,
                                            ProcessInfo& /*rCurrentProcessInfo*/)
{
    const SizeType NumNodes(3), LocalSize(9);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    IndexType LocalIndex = 0;

    for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_SCALAR_1).EquationId();
    }
}

template <>
void VMSAdjointElement<3>::EquationIdVector(EquationIdVectorType& rResult,
                                            ProcessInfo& /*rCurrentProcessInfo*/)
{
    const SizeType NumNodes(4), LocalSize(16);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    IndexType LocalIndex = 0;

    for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_Y).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_Z).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_SCALAR_1).EquationId();
    }
}

template <>
double VMSAdjointElement<2>::CalculateElementSize(const double Area)
{
    // Diameter of a circle of given Area.
    return 1.128379167 * std::sqrt(Area);
}

template <>
double VMSAdjointElement<3>::CalculateElementSize(const double Volume)
{
    // Diameter of the sphere circumscribed to a regular tetrahedron
    // with the same volume.
    return 0.60046878 * std::pow(Volume, 0.333333333333333333333);
}

template <>
void VMSAdjointElement<2>::CalculateDeterminantOfJacobianDerivatives(array_1d<double, 6>& rDetJDerivatives)
{
    const SizeType LocalSize(6);

    if (rDetJDerivatives.size() != LocalSize)
        rDetJDerivatives.resize(LocalSize, false);

    rDetJDerivatives[0] = this->GetGeometry()[1].Y() - this->GetGeometry()[2].Y();
    rDetJDerivatives[1] = this->GetGeometry()[2].X() - this->GetGeometry()[1].X();
    rDetJDerivatives[2] = this->GetGeometry()[2].Y() - this->GetGeometry()[0].Y();
    rDetJDerivatives[3] = this->GetGeometry()[0].X() - this->GetGeometry()[2].X();
    rDetJDerivatives[4] = this->GetGeometry()[0].Y() - this->GetGeometry()[1].Y();
    rDetJDerivatives[5] = this->GetGeometry()[1].X() - this->GetGeometry()[0].X();
}

template <>
void VMSAdjointElement<3>::CalculateDeterminantOfJacobianDerivatives(array_1d<double, 12>& rDetJDerivatives)
{
    const SizeType LocalSize(12);

    if (rDetJDerivatives.size() != LocalSize)
        rDetJDerivatives.resize(LocalSize, false);

    const double x0 = this->GetGeometry()[0].X();
    const double x1 = this->GetGeometry()[1].X();
    const double x2 = this->GetGeometry()[2].X();
    const double x3 = this->GetGeometry()[3].X();

    const double y0 = this->GetGeometry()[0].Y();
    const double y1 = this->GetGeometry()[1].Y();
    const double y2 = this->GetGeometry()[2].Y();
    const double y3 = this->GetGeometry()[3].Y();

    const double z0 = this->GetGeometry()[0].Z();
    const double z1 = this->GetGeometry()[1].Z();
    const double z2 = this->GetGeometry()[2].Z();
    const double z3 = this->GetGeometry()[3].Z();

    rDetJDerivatives[0] = -z1 * y3 + z1 * y2 + y3 * z2 - y2 * z3 + y1 * z3 - y1 * z2;
    rDetJDerivatives[1] = -z1 * x2 + z1 * x3 + x1 * z2 - z2 * x3 + x2 * z3 - x1 * z3;
    rDetJDerivatives[2] = -x2 * y3 + y2 * x3 - y1 * x3 + x1 * y3 + y1 * x2 - x1 * y2;
    rDetJDerivatives[3] = y2 * z3 - y2 * z0 - y0 * z3 - y3 * z2 + y3 * z0 + y0 * z2;
    rDetJDerivatives[4] = z2 * x3 - z2 * x0 - z0 * x3 - x2 * z3 + x2 * z0 + x0 * z3;
    rDetJDerivatives[5] = x2 * y3 - x2 * y0 - x0 * y3 - y2 * x3 + y2 * x0 + y0 * x3;
    rDetJDerivatives[6] = z1 * y3 - z1 * y0 - y3 * z0 - y1 * z3 + y0 * z3 + y1 * z0;
    rDetJDerivatives[7] = -z1 * x3 + z1 * x0 + z0 * x3 - x0 * z3 + x1 * z3 - x1 * z0;
    rDetJDerivatives[8] = x1 * y0 + x0 * y3 - y0 * x3 - x1 * y3 + y1 * x3 - y1 * x0;
    rDetJDerivatives[9] = -z1 * y2 + z1 * y0 + y2 * z0 - y0 * z2 + y1 * z2 - y1 * z0;
    rDetJDerivatives[10] = z1 * x2 - z1 * x0 - x2 * z0 + z2 * x0 - x1 * z2 + x1 * z0;
    rDetJDerivatives[11] = -y2 * x0 + x1 * y2 + y1 * x0 - y1 * x2 + x2 * y0 - x1 * y0;
}

template <>
void VMSAdjointElement<2>::CalculateStabilizationParametersDerivative(
    double& rTauOneDeriv,
    double& rTauTwoDeriv,
    const double TauOne,
    const double TauTwo,
    const double VelNorm,
    const double ElemSize,
    const double Density,
    const double Viscosity,
    const double DetJDeriv)
{
    const double TwoPi = 6.283185307179586;
    const double TwoOverPi = 0.636619772367581;
    rTauOneDeriv = TwoOverPi * TauOne * TauOne / std::pow(ElemSize, 3) *
                   (Density * VelNorm + 4.0 * Viscosity / ElemSize) * DetJDeriv;

    rTauTwoDeriv = Density * VelNorm / (TwoPi * ElemSize) * DetJDeriv;
}

template <>
void VMSAdjointElement<3>::CalculateStabilizationParametersDerivative(
    double& rTauOneDeriv,
    double& rTauTwoDeriv,
    const double TauOne,
    const double TauTwo,
    const double VelNorm,
    const double ElemSize,
    const double Density,
    const double Viscosity,
    const double DetJDeriv)
{
    const double CoefOne = 0.0240562975623840;
    const double CoefTwo = 0.00601407439059599;

    rTauOneDeriv = CoefOne * TauOne * TauOne / std::pow(ElemSize, 4) *
                   (Density * VelNorm + 4.0 * Viscosity / ElemSize) * DetJDeriv;

    rTauTwoDeriv = CoefTwo * Density * VelNorm / (ElemSize * ElemSize) * DetJDeriv;
}

template <>
void VMSAdjointElement<2>::AddViscousTerm(MatrixType& rResult,
                                          const VMSAdjointElement<2>::ShapeFunctionDerivativesType& rDN_DX,
                                          const double Weight)
{
    const SizeType NumNodes = 3;

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    IndexType FirstRow(0), FirstCol(0);

    for (IndexType j = 0; j < NumNodes; ++j)
    {
        for (IndexType i = 0; i < NumNodes; ++i)
        {
            // First Row
            rResult(FirstRow, FirstCol) +=
                Weight * (FourThirds * rDN_DX(i, 0) * rDN_DX(j, 0) +
                          rDN_DX(i, 1) * rDN_DX(j, 1));
            rResult(FirstRow, FirstCol + 1) +=
                Weight * (nTwoThirds * rDN_DX(i, 0) * rDN_DX(j, 1) +
                          rDN_DX(i, 1) * rDN_DX(j, 0));

            // Second Row
            rResult(FirstRow + 1, FirstCol) +=
                Weight * (nTwoThirds * rDN_DX(i, 1) * rDN_DX(j, 0) +
                          rDN_DX(i, 0) * rDN_DX(j, 1));
            rResult(FirstRow + 1, FirstCol + 1) +=
                Weight * (FourThirds * rDN_DX(i, 1) * rDN_DX(j, 1) +
                          rDN_DX(i, 0) * rDN_DX(j, 0));

            // Update Counter
            FirstRow += 3;
        }
        FirstRow = 0;
        FirstCol += 3;
    }
}

template <>
void VMSAdjointElement<3>::AddViscousTerm(MatrixType& rResult,
                                          const VMSAdjointElement<3>::ShapeFunctionDerivativesType& rDN_DX,
                                          const double Weight)
{
    const unsigned int NumNodes = 4;

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0), FirstCol(0);

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag = rDN_DX(i, 0) * rDN_DX(j, 0) +
                                rDN_DX(i, 1) * rDN_DX(j, 1) +
                                rDN_DX(i, 2) * rDN_DX(j, 2);

            // First Row
            rResult(FirstRow, FirstCol) +=
                Weight * (OneThird * rDN_DX(i, 0) * rDN_DX(j, 0) + Diag);
            rResult(FirstRow, FirstCol + 1) +=
                Weight * (nTwoThirds * rDN_DX(i, 0) * rDN_DX(j, 1) +
                          rDN_DX(i, 1) * rDN_DX(j, 0));
            rResult(FirstRow, FirstCol + 2) +=
                Weight * (nTwoThirds * rDN_DX(i, 0) * rDN_DX(j, 2) +
                          rDN_DX(i, 2) * rDN_DX(j, 0));

            // Second Row
            rResult(FirstRow + 1, FirstCol) +=
                Weight * (nTwoThirds * rDN_DX(i, 1) * rDN_DX(j, 0) +
                          rDN_DX(i, 0) * rDN_DX(j, 1));
            rResult(FirstRow + 1, FirstCol + 1) +=
                Weight * (OneThird * rDN_DX(i, 1) * rDN_DX(j, 1) + Diag);
            rResult(FirstRow + 1, FirstCol + 2) +=
                Weight * (nTwoThirds * rDN_DX(i, 1) * rDN_DX(j, 2) +
                          rDN_DX(i, 2) * rDN_DX(j, 1));

            // Third Row
            rResult(FirstRow + 2, FirstCol) +=
                Weight * (nTwoThirds * rDN_DX(i, 2) * rDN_DX(j, 0) +
                          rDN_DX(i, 0) * rDN_DX(j, 2));
            rResult(FirstRow + 2, FirstCol + 1) +=
                Weight * (nTwoThirds * rDN_DX(i, 2) * rDN_DX(j, 1) +
                          rDN_DX(i, 1) * rDN_DX(j, 2));
            rResult(FirstRow + 2, FirstCol + 2) +=
                Weight * (OneThird * rDN_DX(i, 2) * rDN_DX(j, 2) + Diag);

            // Update Counter
            FirstRow += 4;
        }
        FirstRow = 0;
        FirstCol += 4;
    }
}

template <>
void VMSAdjointElement<2>::AddViscousTermDerivative(
    BoundedMatrix<double, 9, 9>& rResult,
    const VMSAdjointElement<2>::ShapeFunctionDerivativesType& rDN_DX,
    const VMSAdjointElement<2>::ShapeFunctionDerivativesType& rDN_DX_Deriv,
    const double Weight,
    const double WeightDeriv)
{
    const SizeType NumNodes = 3;

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    IndexType FirstRow(0), FirstCol(0);

    for (IndexType j = 0; j < NumNodes; ++j)
    {
        for (IndexType i = 0; i < NumNodes; ++i)
        {
            // First Row
            rResult(FirstRow, FirstCol) +=
                WeightDeriv * (FourThirds * rDN_DX(i, 0) * rDN_DX(j, 0) +
                               rDN_DX(i, 1) * rDN_DX(j, 1)) +
                Weight * (FourThirds * rDN_DX_Deriv(i, 0) * rDN_DX(j, 0) +
                          rDN_DX_Deriv(i, 1) * rDN_DX(j, 1)) +
                Weight * (FourThirds * rDN_DX(i, 0) * rDN_DX_Deriv(j, 0) +
                          rDN_DX(i, 1) * rDN_DX_Deriv(j, 1));
            rResult(FirstRow, FirstCol + 1) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 0) * rDN_DX(j, 1) +
                               rDN_DX(i, 1) * rDN_DX(j, 0)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 0) * rDN_DX(j, 1) +
                          rDN_DX_Deriv(i, 1) * rDN_DX(j, 0)) +
                Weight * (nTwoThirds * rDN_DX(i, 0) * rDN_DX_Deriv(j, 1) +
                          rDN_DX(i, 1) * rDN_DX_Deriv(j, 0));

            // Second Row
            rResult(FirstRow + 1, FirstCol) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 1) * rDN_DX(j, 0) +
                               rDN_DX(i, 0) * rDN_DX(j, 1)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 1) * rDN_DX(j, 0) +
                          rDN_DX_Deriv(i, 0) * rDN_DX(j, 1)) +
                Weight * (nTwoThirds * rDN_DX(i, 1) * rDN_DX_Deriv(j, 0) +
                          rDN_DX(i, 0) * rDN_DX_Deriv(j, 1));
            rResult(FirstRow + 1, FirstCol + 1) +=
                WeightDeriv * (FourThirds * rDN_DX(i, 1) * rDN_DX(j, 1) +
                               rDN_DX(i, 0) * rDN_DX(j, 0)) +
                Weight * (FourThirds * rDN_DX_Deriv(i, 1) * rDN_DX(j, 1) +
                          rDN_DX_Deriv(i, 0) * rDN_DX(j, 0)) +
                Weight * (FourThirds * rDN_DX(i, 1) * rDN_DX_Deriv(j, 1) +
                          rDN_DX(i, 0) * rDN_DX_Deriv(j, 0));

            // Update Counter
            FirstRow += 3;
        }
        FirstRow = 0;
        FirstCol += 3;
    }
}

template <>
void VMSAdjointElement<3>::AddViscousTermDerivative(
    BoundedMatrix<double, 16, 16>& rResult,
    const VMSAdjointElement<3>::ShapeFunctionDerivativesType& rDN_DX,
    const VMSAdjointElement<3>::ShapeFunctionDerivativesType& rDN_DX_Deriv,
    const double Weight,
    const double WeightDeriv)
{
    const unsigned int NumNodes = 4;

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0), FirstCol(0);

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag = rDN_DX(i, 0) * rDN_DX(j, 0) +
                                rDN_DX(i, 1) * rDN_DX(j, 1) +
                                rDN_DX(i, 2) * rDN_DX(j, 2);
            const double DiagDerivNi = rDN_DX_Deriv(i, 0) * rDN_DX(j, 0) +
                                       rDN_DX_Deriv(i, 1) * rDN_DX(j, 1) +
                                       rDN_DX_Deriv(i, 2) * rDN_DX(j, 2);
            const double DiagDerivNj = rDN_DX(i, 0) * rDN_DX_Deriv(j, 0) +
                                       rDN_DX(i, 1) * rDN_DX_Deriv(j, 1) +
                                       rDN_DX(i, 2) * rDN_DX_Deriv(j, 2);

            // First Row
            rResult(FirstRow, FirstCol) +=
                WeightDeriv * (OneThird * rDN_DX(i, 0) * rDN_DX(j, 0) + Diag) +
                Weight * (OneThird * rDN_DX_Deriv(i, 0) * rDN_DX(j, 0) + DiagDerivNi) +
                Weight * (OneThird * rDN_DX(i, 0) * rDN_DX_Deriv(j, 0) + DiagDerivNj);
            rResult(FirstRow, FirstCol + 1) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 0) * rDN_DX(j, 1) +
                               rDN_DX(i, 1) * rDN_DX(j, 0)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 0) * rDN_DX(j, 1) +
                          rDN_DX_Deriv(i, 1) * rDN_DX(j, 0)) +
                Weight * (nTwoThirds * rDN_DX(i, 0) * rDN_DX_Deriv(j, 1) +
                          rDN_DX(i, 1) * rDN_DX_Deriv(j, 0));
            rResult(FirstRow, FirstCol + 2) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 0) * rDN_DX(j, 2) +
                               rDN_DX(i, 2) * rDN_DX(j, 0)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 0) * rDN_DX(j, 2) +
                          rDN_DX_Deriv(i, 2) * rDN_DX(j, 0)) +
                Weight * (nTwoThirds * rDN_DX(i, 0) * rDN_DX_Deriv(j, 2) +
                          rDN_DX(i, 2) * rDN_DX_Deriv(j, 0));

            // Second Row
            rResult(FirstRow + 1, FirstCol) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 1) * rDN_DX(j, 0) +
                               rDN_DX(i, 0) * rDN_DX(j, 1)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 1) * rDN_DX(j, 0) +
                          rDN_DX_Deriv(i, 0) * rDN_DX(j, 1)) +
                Weight * (nTwoThirds * rDN_DX(i, 1) * rDN_DX_Deriv(j, 0) +
                          rDN_DX(i, 0) * rDN_DX_Deriv(j, 1));
            rResult(FirstRow + 1, FirstCol + 1) +=
                WeightDeriv * (OneThird * rDN_DX(i, 1) * rDN_DX(j, 1) + Diag) +
                Weight * (OneThird * rDN_DX_Deriv(i, 1) * rDN_DX(j, 1) + DiagDerivNi) +
                Weight * (OneThird * rDN_DX(i, 1) * rDN_DX_Deriv(j, 1) + DiagDerivNj);
            rResult(FirstRow + 1, FirstCol + 2) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 1) * rDN_DX(j, 2) +
                               rDN_DX(i, 2) * rDN_DX(j, 1)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 1) * rDN_DX(j, 2) +
                          rDN_DX_Deriv(i, 2) * rDN_DX(j, 1)) +
                Weight * (nTwoThirds * rDN_DX(i, 1) * rDN_DX_Deriv(j, 2) +
                          rDN_DX(i, 2) * rDN_DX_Deriv(j, 1));

            // Third Row
            rResult(FirstRow + 2, FirstCol) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 2) * rDN_DX(j, 0) +
                               rDN_DX(i, 0) * rDN_DX(j, 2)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 2) * rDN_DX(j, 0) +
                          rDN_DX_Deriv(i, 0) * rDN_DX(j, 2)) +
                Weight * (nTwoThirds * rDN_DX(i, 2) * rDN_DX_Deriv(j, 0) +
                          rDN_DX(i, 0) * rDN_DX_Deriv(j, 2));
            rResult(FirstRow + 2, FirstCol + 1) +=
                WeightDeriv * (nTwoThirds * rDN_DX(i, 2) * rDN_DX(j, 1) +
                               rDN_DX(i, 1) * rDN_DX(j, 2)) +
                Weight * (nTwoThirds * rDN_DX_Deriv(i, 2) * rDN_DX(j, 1) +
                          rDN_DX_Deriv(i, 1) * rDN_DX(j, 2)) +
                Weight * (nTwoThirds * rDN_DX(i, 2) * rDN_DX_Deriv(j, 1) +
                          rDN_DX(i, 1) * rDN_DX_Deriv(j, 2));
            rResult(FirstRow + 2, FirstCol + 2) +=
                WeightDeriv * (OneThird * rDN_DX(i, 2) * rDN_DX(j, 2) + Diag) +
                Weight * (OneThird * rDN_DX_Deriv(i, 2) * rDN_DX(j, 2) + DiagDerivNi) +
                Weight * (OneThird * rDN_DX(i, 2) * rDN_DX_Deriv(j, 2) + DiagDerivNj);

            // Update Counter
            FirstRow += 4;
        }
        FirstRow = 0;
        FirstCol += 4;
    }
}

///@} // Specialized implementations

template class VMSAdjointElement<2>;
template class VMSAdjointElement<3>;

} // namespace Kratos
