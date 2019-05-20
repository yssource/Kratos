//  KratosFluidDynamicsApplication
//
//  License:		 BSD License
//					 license: FluidDynamicsApplication/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/indirect_scalar.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief An adjoint element for discrete shape sensitivity of VMS fluid element.
 *
 * @see VMS monolithic fluid element
 */
template <unsigned int TDim>
class VMSAdjointElement : public Element
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement);

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step) override;

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step) override;

        void GetAuxiliaryVector(std::size_t NodeId,
                                std::vector<IndirectScalar<double>>& rVector,
                                std::size_t Step) override;

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override;
    };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(VMSAdjointElement);

    constexpr static unsigned int TNumNodes = TDim + 1;

    constexpr static unsigned int TBlockSize = TDim + 1;

    constexpr static unsigned int TFluidLocalSize = TBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    typedef Element::IndexType IndexType;

    typedef Element::SizeType SizeType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::VectorType VectorType;

    typedef Element::MatrixType MatrixType;

    typedef Element::DofsVectorType DofsVectorType;

    typedef Element::EquationIdVectorType EquationIdVectorType;

    typedef BoundedMatrix<double, TNumNodes, TDim> ShapeFunctionDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    VMSAdjointElement(IndexType NewId = 0);

    VMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry);

    VMSAdjointElement(IndexType NewId,
                      GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);

    ~VMSAdjointElement() override;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    /**
     * @brief Creates a new element of this type.
     *
     * @return pointer to the newly created element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Checks for proper element geometry, nodal variables and dofs.
     *
     * @return 0 after successful completion.
     */
    int Check(const ProcessInfo& /*rCurrentProcessInfo*/) override;

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(VectorType& rValues, int Step = 0) override;

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override;

    /// Returns the adjoint acceleration values stored in this element's nodes.
    void GetSecondDerivativesVector(VectorType& rValues, int Step = 0) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& /*rCurrentProcessInfo*/) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& /*rCurrentProcessInfo*/) override;

    /**
     * @brief Calculates the adjoint matrix for velocity and pressure.
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * velocity and pressure transposed:
     *
     * \f[
     *    \partial_{\mathbf{w}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{w}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     *
     * where \f$\mathbf{w}^n\f$ is the vector of nodal velocities and pressures
     * stored at the current step. For steady problems, the ACCELERATION
     * (\f$\dot{\mathbf{w}}^n\f$) must be set to zero on the nodes. For
     * the Bossak method, \f$\dot{\mathbf{w}}^{n-\alpha}\f$ must be stored in
     * ACCELERATION.
     */
    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculates the adjoint matrix for acceleration.
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * acceleration:
     *
     * \f[
     *    \partial_{\dot{\mathbf{w}}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\dot{\mathbf{w}}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             ProcessInfo& /*rCurrentProcessInfo*/) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& /*rCurrentProcessInfo*/) override;

    /**
     * @brief Calculates the sensitivity matrix.
     *
     * \f[
     *    \partial_{\mathbf{s}}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{s}}(\mathbf{M}^n \dot{\mathbf{w}}^{n-\alpha})^T
     * \f]
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    ProcessInfo& /*rCurrentProcessInfo*/) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& /*rCurrentProcessInfo*/) override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    /// Calculate VMS-stabilized (lumped) mass matrix.
    void CalculateVMSMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Adds primal gradient of the VMS mass matrix multiplied by a vector.
     *
     * Calculates \f$ \partial_{\mathbf{w}^n} (\mathbf{M}^n\mathbf{x}) \f$.
     * \f$\mathbf{w}^n\f$ is the vector of primal variables at the current
     * adjoint step. \f$\mathbf{M}^n\f$ is the VMS mass matrix. \f$\mathbf{x}\f$
     * is a constant vector with zeros for pressure dofs. The variable
     * determines the values for velocity dofs.
     */
    void AddPrimalGradientOfVMSMassTerm(MatrixType& rOutputMatrix,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        double alpha,
                                        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Adds shape gradient of the VMS mass matrix multiplied by a vector.
     *
     * Calculates \f$ \partial_{\mathbf{s}} (\mathbf{M}^n\mathbf{x})^T \f$.
     * \f$\mathbf{s}\f$ is the vector of nodal coordinates. \f$\mathbf{M}^n\f$.
     * is the VMS mass matrix at the current adjoint step. \f$\mathbf{x}\f$ is
     * a constant vector with zeros for pressure dofs. The variable
     * determines the values for velocity dofs.
     */
    void AddShapeGradientOfVMSMassTerm(MatrixType& rOutputMatrix,
                                       const Variable<array_1d<double, 3>>& rVariable,
                                       double alpha,
                                       const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculates the elemental contribution to the steady adjoint system matrix.
     *
     * This function returns elemental contributions for:
     *
     * \f[
     * \partial_{\mathbf{w}^n}\mathbf{f}(\mathbf{w}^n)
     * \f]
     *
     * where the current adjoint step is the \f$n^{th}\f$ time step.
     */
    void CalculatePrimalGradientOfVMSSteadyTerm(MatrixType& rAdjointMatrix,
                                                const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate the partial derivatives of damping term w.r.t. shape parameters.
     *
     * This function returns elemental contributions for:
     *
     * \f[
     * \partial_{\mathbf{s}}\mathbf{f}(\mathbf{w}^n)^T
     * \f]
     *
     * \f$\mathbf{s}\f$ are the coordinates of the element's nodes.
     *
     * This function is only valid when the determinant of the Jacobian is constant
     * over the element.
     */
    void CalculateShapeGradientOfVMSSteadyTerm(MatrixType& rShapeDerivativesMatrix,
                                               const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Returns the gradient matrix of the velocity.
     *
     * The row index corresponds to the velocity component and the column index
     * to the derivative.
     *
     * @param rGradVel velocity gradient matrix
     * @param rDN_DX shape functions' gradients
     */
    void CalculateVelocityGradient(BoundedMatrix<double, TDim, TDim>& rGradVel,
                                   const ShapeFunctionDerivativesType& rDN_DX);

    /**
     * @brief Returns the pressure gradient.
     *
     * @param rGradP pressure gradient
     * @param rDN_DX shape functions' gradients
     */
    void CalculatePressureGradient(array_1d<double, TDim>& rGradP,
                                   const ShapeFunctionDerivativesType& rDN_DX);

    /**
     * @brief Returns the element's size.
     *
     * @param Volume the volume (area in 2D) of the element
     */
    double CalculateElementSize(const double Volume);

    /**
     * @brief Returns derivatives of determinant of Jacobian w.r.t coordinates.
     *
     * The derivative of the determinant of the Jacobian w.r.t the jth
     * coordinate of the ith node is stored at the index (i * TDim + j).
     *
     * This function is only valid when the determinant of the Jacobian is
     * constant over the element.
     *
     * @see Triangle2D3
     * @see Tetrahedra3D4
     */
    void CalculateDeterminantOfJacobianDerivatives(array_1d<double, TCoordLocalSize>& rDetJDerivatives);

    /**
     * @brief Returns the VMS stabilization parameters.
     *
     * @param rTauOne momentum stabilization parameter
     * @param rTauTwo divergence stabilization parameter
     * @param VelNorm Euclidean norm of the velocity
     * @param ElemSize size of this element
     * @param Density density of the fluid
     * @param Viscosity dynamic viscosity of the fluid
     */
    void CalculateStabilizationParameters(double& rTauOne,
                                          double& rTauTwo,
                                          double VelNorm,
                                          double ElemSize,
                                          double Density,
                                          double Viscosity,
                                          const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Returns stabilization parameters derived w.r.t a node's coordinate.
     *
     * @param rTauOneDeriv derivative of momentum stabilization parameter
     * @param rTauTwoDeriv derivative of divergence stabilization parameter
     * @param TauOne momentum stabilization parameter
     * @param TauTwo divergence stabilization parameter
     * @param VelNorm Euclidean norm of the velocity
     * @param ElemSize size of this element
     * @param Density density of the fluid
     * @param Viscosity dynamic viscosity of the fluid
     * @param DetJDeriv derivative of the determinant of the Jacobian
     */
    void CalculateStabilizationParametersDerivative(double& rTauOneDeriv,
                                                    double& rTauTwoDeriv,
                                                    double TauOne,
                                                    double TauTwo,
                                                    double VelNorm,
                                                    double ElemSize,
                                                    double Density,
                                                    double Viscosity,
                                                    double DetJDeriv);

    /**
     * @brief Returns a scalar variable at this integration point.
     *
     * @param rResult the value of the scalar variable at this integration point
     * @param rVariable the variable to be evaluated
     * @param rShapeFunc array of shape function values at this integration point
     */
    void EvaluateInPoint(double& rResult,
                         const Variable<double>& rVariable,
                         const array_1d<double, TNumNodes>& rShapeFunc,
                         IndexType step = 0);

    /**
     * @brief Returns a vector variable at this integration point.
     *
     * @param rResult the value of the vector variable at this integration point
     * @param rVariable the variable to be evaluated
     * @param rShapeFunc array of shape function values at this integration point
     */
    void EvaluateInPoint(array_1d<double, TDim>& rResult,
                         const Variable<array_1d<double, 3>>& rVariable,
                         const array_1d<double, TNumNodes>& rN,
                         IndexType step = 0);

    /**
     * @brief Adds viscous contributions to adjoint system matrix.
     *
     * @param rResult matrix to add viscous contributions to
     * @param rDN_DX shape functions' gradients
     * @param Weight integration weight including dynamic viscosity
     */
    void AddViscousTerm(MatrixType& rResult,
                        const ShapeFunctionDerivativesType& rDN_DX,
                        const double Weight);

    /**
     * @brief Adds derivative of viscous term w.r.t a node's coordinate.
     *
     * @param rResult matrix to add viscous contributions to
     * @param rDN_DX shape functions' gradients
     * @param rDN_DX_Deriv shape functions' gradients derived w.r.t the coordinate
     * @param Weight integration weight including dynamic viscosity
     * @param WeightDeriv integration weight derived w.r.t the coordinate
     *
     * @see AddViscousTerm
     */
    void AddViscousTermDerivative(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rResult,
                                  const ShapeFunctionDerivativesType& rDN_DX,
                                  const ShapeFunctionDerivativesType& rDN_DX_Deriv,
                                  const double Weight,
                                  const double WeightDeriv);

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY;

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY;

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Unaccessible methods
    ///@{

    VMSAdjointElement& operator=(VMSAdjointElement const& rOther);

    VMSAdjointElement(VMSAdjointElement const& rOther);

    ///@}

}; // class VMSAdjointElement

///@} // Kratos classes

///@name Input and output
///@{

/// Defines an input stream operator that does nothing.
template <unsigned int TDim>
inline std::istream& operator>>(std::istream& rIStream, VMSAdjointElement<TDim>& rThis)
{
    return rIStream;
}

/// Defines an output stream operator that prints element info.
template <unsigned int TDim>
inline std::ostream& operator<<(std::ostream& rOStream, const VMSAdjointElement<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // FluidDynamicsApplication group

} // namespace Kratos

#endif // KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED defined
