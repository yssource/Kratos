//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

// System includes

// External includes

// Include Base h
#include "custom_elements/rans_constitutive_element.h"
#include "custom_utilities/calculation_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Constructor.
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::RANSConstitutiveElement(IndexType NewId)
    : Element(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::RANSConstitutiveElement(
    IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::RANSConstitutiveElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::RANSConstitutiveElement(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::RANSConstitutiveElement(RANSConstitutiveElement const& rOther)
    : Element(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::~RANSConstitutiveElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>& RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::operator=(
    RANSConstitutiveElement<TDim, TNumNodes, TBlockSize> const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator=(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

///@}
///@name Operations
///@{

/**
 * ELEMENTS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
Element::Pointer RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY;
    KRATOS_ERROR
        << "Attempting to Create base RANSConstitutiveElement instances." << std::endl;
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
Element::Pointer RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY;
    KRATOS_ERROR
        << "Attempting to Create base RANSConstitutiveElement instances." << std::endl;
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
Element::Pointer RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_shared<RANSConstitutiveElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Attempting to call base RANSConstitutiveElement "
                    "EquationIdVector method."
                 << std::endl;
    KRATOS_CATCH("");
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::GetDofList(
    DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR
        << "Attempting to call base RANSConstitutiveElement GetDofList method."
        << std::endl;
    KRATOS_CATCH("");
}

/**
 * ELEMENTS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
 * methods they can be managed internally with a private method to do the same
 * calculations only once: MANDATORY
 */

/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != TLocalSize || rLeftHandSideMatrix.size2() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);

    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TLocalSize, TLocalSize);
    noalias(rRightHandSideVector) = ZeroVector(TLocalSize);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != TLocalSize || rLeftHandSideMatrix.size2() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TLocalSize, TLocalSize);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(TLocalSize);
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rDampMatrix.size1() != TLocalSize)
        rDampMatrix.resize(TLocalSize, TLocalSize, false);

    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);

    noalias(rDampMatrix) = ZeroMatrix(TLocalSize, TLocalSize);
    noalias(rRightHandSideVector) = ZeroVector(TLocalSize);
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TLocalSize || rLeftHandSideMatrix.size2() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);
    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TLocalSize || rLeftHandSideMatrix.size2() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateFirstDerivativesRHS(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);
}

/**
 * ELEMENTS inherited from this class must implement this methods
 * if they need to add dynamic element contributions
 * note: second derivatives means the accelerations if the displacements are the dof of the analysis
 * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
 * CalculateSecondDerivativesContributions,
 * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
 */

/**
 * this is called during the assembling process in order
 * to calculate the second derivative contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0 || rLeftHandSideMatrix.size2() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);
    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TLocalSize || rLeftHandSideMatrix.size2() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateMassMatrix(
    MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != TLocalSize || rMassMatrix.size2() != TLocalSize)
        rMassMatrix.resize(TLocalSize, TLocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(TLocalSize, TLocalSize);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != TLocalSize || rDampingMatrix.size2() != TLocalSize)
        rDampingMatrix.resize(TLocalSize, TLocalSize, false);

    noalias(rDampingMatrix) = ZeroMatrix(TLocalSize, TLocalSize);
}

/**
 * This method provides the place to perform checks on the completeness of the
 * input and the compatibility with the problem options as well as the
 * contitutive laws selected It is designed to be called only once (or anyway,
 * not often) typically at the beginning of the calculations, so to verify that
 * nothing is missing from the input or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
int RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)
        << "RANSConstitutiveElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0)
        << "On RANSConstitutiveElement -> " << this->Id()
        << "; Area cannot be less than or equal to 0" << std::endl;

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
GeometryData::IntegrationMethod RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateGeometryData(
    Vector& rGaussWeights, Matrix& rNContainer, ShapeFunctionDerivativesArrayType& rDN_DX) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    CalculationUtilities::CalculateGeometryData(
        r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
double RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::EvaluateInPoint(
    const Variable<double>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return CalculationUtilities::EvaluateInPoint<Geometry<Node<3>>>(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
array_1d<double, 3> RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::EvaluateInPoint(
    const Variable<array_1d<double, 3>>& rVariable, const Vector& rShapeFunction, const int Step) const
{
    return CalculationUtilities::EvaluateInPoint<Geometry<Node<3>>>(
        this->GetGeometry(), rVariable, rShapeFunction, Step);
}
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::GetConvectionOperator(
    BoundedVector<double, TNumNodes>& rOutput,
    const array_1d<double, 3>& rVector,
    const Matrix& rShapeDerivatives) const
{
    rOutput.clear();
    for (unsigned int i = 0; i < TNumNodes; ++i)
        for (unsigned int j = 0; j < TDim; j++)
        {
            rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
        }
}
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
double RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::GetDivergenceOperator(
    const Variable<array_1d<double, 3>>& rVariable, const Matrix& rShapeDerivatives, const int Step) const
{
    double value = 0.0;
    const GeometryType& r_geometry = this->GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const array_1d<double, 3>& r_value =
            r_geometry[i].FastGetSolutionStepValue(rVariable, Step);
        for (unsigned int j = 0; j < TDim; ++j)
        {
            value += r_value[j] * rShapeDerivatives(i, j);
        }
    }

    return value;
}
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateGradient(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step) const
{
    rOutput.clear();
    const GeometryType& r_geometry = this->GetGeometry();

    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        const array_1d<double, 3>& r_value =
            r_geometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (unsigned int i = 0; i < TDim; ++i)
        {
            for (unsigned int j = 0; j < TDim; ++j)
            {
                rOutput(i, j) += rShapeDerivatives(a, j) * r_value[i];
            }
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateGradient(
    array_1d<double, 3>& rOutput,
    const Variable<double>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step) const
{
    rOutput.clear();
    const GeometryType& r_geometry = this->GetGeometry();
    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        const double value = r_geometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (unsigned int i = 0; i < TDim; ++i)
            rOutput[i] += rShapeDerivatives(a, i) * value;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateSymmetricGradientMatrix(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rVariable,
    const BoundedMatrix<double, TDim, TDim>& rGradientMatrix,
    const Matrix& rShapeDerivatives,
    const int Step) const
{
    const double variable_divergence =
        this->GetDivergenceOperator(rVariable, rShapeDerivatives, Step);
    identity_matrix<double> identity(TDim);

    rOutput.clear();
    noalias(rOutput) = rGradientMatrix + trans(rGradientMatrix) -
                       (2.0 / 3.0) * variable_divergence * identity;
}

///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a string.
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
std::string RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::Info() const
{
    std::stringstream buffer;
    buffer << "RANSConstitutiveElement #" << Id();
    return buffer.str();
}

/// Print information about this object.
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RANSConstitutiveElement #" << Id();
}

/// Print object's data.
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

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

// template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
// double RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateUTau(
//     const double velocity_magnitude,
//     const double wall_distance,
//     const double nu,
//     const double beta,
//     const double von_karman) const
// {
//     // auto log_law_wall_function = [velocity_magnitude, von_karman,
//     //                               wall_distance, nu, beta](double u_tau) {
//     //     return u_tau -
//     //            velocity_magnitude /
//     //                ((1.0 / von_karman) * std::log(u_tau * wall_distance / nu) + beta);
//     // };

//     // const double iterTol = 1e-5;
//     // const double maxIter = 10;
//     // u_tau = BrentIteration::FindRoot(log_law_wall_function, 1e-6,
//     //                                  10.0 * nu / wall_distance, iterTol, maxIter);

//     const unsigned int max_iterations = 10;
//     double u_tau = nu / wall_distance;
//     for (unsigned int i = 0; i < max_iterations; ++i)
//     {
//         u_tau = nu * std::exp(((velocity_magnitude / u_tau) - beta) * von_karman) / wall_distance;
//     }

//     return u_tau;
// }

// template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
// double RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::CalculateYPlus(
//     const double u_tau, const double velocity_magnitude, const double wall_distance, const double nu) const
// {
//     double y_plus = u_tau * wall_distance / nu;
//     if (y_plus < 11.06)
//         y_plus = velocity_magnitude / u_tau;

//     return y_plus;
// }

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

///@name Static Member Variables
///@{

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
///@name Serialization
///@{
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

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
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
inline std::istream& operator>>(std::istream& rIStream,
                                RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RANSConstitutiveElement<TDim, TNumNodes, TBlockSize>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
// Class template instantiation

template class RANSConstitutiveElement<2, 3, 1>;
template class RANSConstitutiveElement<3, 4, 1>;
template class RANSConstitutiveElement<2, 4, 1>;
template class RANSConstitutiveElement<3, 8, 1>;

} // namespace Kratos.
