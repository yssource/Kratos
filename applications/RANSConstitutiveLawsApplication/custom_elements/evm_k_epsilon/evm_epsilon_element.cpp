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
#include "evm_epsilon_element.h"

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
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>::EvmEpsilonElement(IndexType NewId)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>::EvmEpsilonElement(IndexType NewId,
                                                      const NodesArrayType& ThisNodes)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>::EvmEpsilonElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>::EvmEpsilonElement(IndexType NewId,
                                                      GeometryType::Pointer pGeometry,
                                                      PropertiesType::Pointer pProperties)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>::EvmEpsilonElement(EvmEpsilonElement<TDim, TNumNodes> const& rOther)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>::~EvmEpsilonElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template <unsigned int TDim, unsigned int TNumNodes>
EvmEpsilonElement<TDim, TNumNodes>& EvmEpsilonElement<TDim, TNumNodes>::operator=(
    EvmEpsilonElement<TDim, TNumNodes> const& rOther)
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
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer EvmEpsilonElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                            NodesArrayType const& ThisNodes,
                                                            PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmEpsilonElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer EvmEpsilonElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                            GeometryType::Pointer pGeom,
                                                            PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmEpsilonElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer EvmEpsilonElement<TDim, TNumNodes>::Clone(IndexType NewId,
                                                           NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmEpsilonElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                          ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != TLocalSize)
        rResult.resize(TLocalSize, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rResult[i] =
            Element::GetGeometry()[i].GetDof(TURBULENT_ENERGY_DISSIPATION_RATE).EquationId();
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TLocalSize)
        rElementalDofList.resize(TLocalSize);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rElementalDofList[i] =
            Element::GetGeometry()[i].pGetDof(TURBULENT_ENERGY_DISSIPATION_RATE);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
            TURBULENT_ENERGY_DISSIPATION_RATE, Step);
    }
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
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                              VectorType& rRightHandSideVector,
                                                              ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                               ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights, ShapeFunctions, ShapeDerivatives);
    const unsigned int num_gauss_points = GaussWeights.size();

    const GeometryType rGeom = this->GetGeometry();

    const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
    const double von_karman = rCurrentProcessInfo[WALL_VON_KARMAN];
    const double C_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double C1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    const double mixing_length = rCurrentProcessInfo[TURBULENT_MIXING_LENGTH];
    const double turbulent_viscosity_fraction =
        rCurrentProcessInfo[TURBULENT_VISCOSITY_FRACTION];

    bool is_near_wall = false;
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        if (rGeom[i].Is(STRUCTURE))
        {
            is_near_wall = true;
            break;
        }
    }

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const auto& r_shape_derivatives = ShapeDerivatives[g];

        double velocity_divergence = 0.0;
        double epsilon = 0.0;
        double tke = 0.0;
        double wall_distance = 0.0;
        double nu = 0.0;
        array_1d<double, 3> velocity = ZeroVector(3);

        for (unsigned int c = 0; c < TNumNodes; c++)
        {
            const array_1d<double, 3>& r_velocity =
                rGeom[c].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int k = 0; k < TDim; k++)
            {
                velocity_divergence += r_shape_derivatives(c, k) * r_velocity[k];
            }

            epsilon += ShapeFunctions(g, c) *
                       rGeom[c].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            tke += ShapeFunctions(g, c) *
                   rGeom[c].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            wall_distance +=
                ShapeFunctions(g, c) * rGeom[c].FastGetSolutionStepValue(DISTANCE);
            nu += ShapeFunctions(g, c) * rGeom[c].FastGetSolutionStepValue(VISCOSITY);
            velocity += ShapeFunctions(g, c) * r_velocity;
        }

        double y_plus(300.0);
        if (is_near_wall)
        {
            const double velocity_magnitude = norm_2(velocity);
            const double u_tau = this->CalculateUTau(
                velocity_magnitude, wall_distance, nu, beta, von_karman);

            y_plus = this->CalculateYPlus(u_tau, velocity_magnitude, wall_distance, nu);
        }

        const double f_mu = 1 - std::exp(-0.0115 * y_plus);
        // calculating limited mixing length
        double temp_tke = C_mu * f_mu * std::pow(tke, 1.5);
        const double limited_mixing_length =
            (temp_tke < epsilon * mixing_length) ? temp_tke / epsilon : mixing_length;

        const double nu_min = nu * turbulent_viscosity_fraction;
        const double nu_t = (nu_min < limited_mixing_length * std::pow(tke, 0.5))
                                ? limited_mixing_length * std::pow(tke, 0.5)
                                : nu_min;
        const double gamma = C_mu * f_mu * tke / nu_t;

        BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
        velocity_gradient_matrix.clear();
        for (unsigned int c = 0; c < TNumNodes; c++)
        {
            const array_1d<double, 3>& r_velocity =
                rGeom[c].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int i = 0; i < TDim; i++)
            {
                for (unsigned int j = 0; j < TDim; j++)
                {
                    velocity_gradient_matrix(i, j) +=
                        r_shape_derivatives(c, j) * r_velocity[i];
                }
            }
        }

        BoundedMatrix<double, TDim, TDim> temp;
        noalias(temp) = velocity_gradient_matrix + trans(velocity_gradient_matrix);
        const double P_k = 0.5 * nu_t * std::pow(norm_frobenius(temp), 2.0);

        for (unsigned int a = 0; a < TNumNodes; a++)
        {
            const double value = GaussWeights[g] * ShapeFunctions(g, a) * P_k;
            // Add epsilon contribution
            rRightHandSideVector[a] += gamma * C1 * value;
        }
    }

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                                      ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
                                                                      ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
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
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                             ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_TRY

    BaseType::CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights, ShapeFunctions, ShapeDerivatives);
    const unsigned int num_gauss_points = GaussWeights.size();

    for (unsigned int a = 0; a < TNumNodes; a++)
    {
        for (unsigned int b = 0; b < TNumNodes; b++)
        {
            for (unsigned int g = 0; g < num_gauss_points; g++)
            {
                const double value =
                    GaussWeights[g] * ShapeFunctions(g, a) * ShapeFunctions(g, b);

                // Add contribution for k
                rMassMatrix(a, b) += value;
            }
        }
    }

    // KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights, ShapeFunctions, ShapeDerivatives);
    const unsigned int num_gauss_points = GaussWeights.size();

    const GeometryType rGeom = this->GetGeometry();

    const double beta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
    const double von_karman = rCurrentProcessInfo[WALL_VON_KARMAN];
    const double C_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double C2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];
    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double mixing_length = rCurrentProcessInfo[TURBULENT_MIXING_LENGTH];
    const double turbulent_viscosity_fraction =
        rCurrentProcessInfo[TURBULENT_VISCOSITY_FRACTION];

    bool is_near_wall = false;
    for (unsigned int i = 0; i < TNumNodes; ++i)
        if (rGeom[i].Is(STRUCTURE))
        {
            is_near_wall = true;
            break;
        }

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const auto& r_shape_derivatives = ShapeDerivatives[g];

        double velocity_divergence = 0.0;
        double epsilon = 0.0;
        double tke = 0.0;
        double wall_distance = 0.0;
        double nu = 0.0;
        array_1d<double, 3> velocity = ZeroVector(3);

        for (unsigned int c = 0; c < TNumNodes; c++)
        {
            const array_1d<double, 3>& r_velocity =
                rGeom[c].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int k = 0; k < TDim; k++)
            {
                velocity_divergence += r_shape_derivatives(c, k) * r_velocity[k];
            }

            epsilon += ShapeFunctions(g, c) *
                       rGeom[c].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            tke += ShapeFunctions(g, c) *
                   rGeom[c].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            wall_distance +=
                ShapeFunctions(g, c) * rGeom[c].FastGetSolutionStepValue(DISTANCE);
            nu += ShapeFunctions(g, c) * rGeom[c].FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            velocity += ShapeFunctions(g, c) * r_velocity;
        }

        double y_plus(300.0);
        if (is_near_wall)
        {
            const double velocity_magnitude = norm_2(velocity);
            const double u_tau = this->CalculateUTau(
                velocity_magnitude, wall_distance, nu, beta, von_karman);
            y_plus = this->CalculateYPlus(u_tau, velocity_magnitude, wall_distance, nu);
        }

        const double f_mu = 1 - std::exp(-0.0115 * y_plus);
        const double f2 =
            1 - 0.22 * std::exp(-std::pow(std::pow(tke, 2) / (6.0 * nu * epsilon), 2));

        const double limited_mixing_length =
            std::min<double>(C_mu * std::pow(tke, 1.5) / epsilon, mixing_length);

        const double nu_min = nu * turbulent_viscosity_fraction;
        const double nu_t =
            std::max<double>(nu_min, limited_mixing_length * std::pow(tke, 0.5));
        const double gamma = C_mu * tke / (f_mu * nu_t);

        for (unsigned int a = 0; a < TNumNodes; a++)
        {
            for (unsigned int b = 0; b < TNumNodes; b++)
            {
                double dNa_dNb = 0.0;
                for (unsigned int i = 0; i < TDim; i++)
                    dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);

                double value_epsilon = 0.0;

                value_epsilon += velocity_divergence;
                value_epsilon += C2 * f2 * gamma;
                value_epsilon += (2.0 * nu / std::pow(wall_distance, 2.0)) *
                                 std::exp(-0.5 * y_plus);
                value_epsilon *= ShapeFunctions(g, a) * ShapeFunctions(g, b);
                value_epsilon += (nu_t / epsilon_sigma) * dNa_dNb;

                // Add epsilon contribution
                rDampingMatrix(a, b) += GaussWeights[g] * value_epsilon;
            }
        }
    }

    KRATOS_CATCH("");
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
template <unsigned int TDim, unsigned int TNumNodes>
int EvmEpsilonElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)
        << "EvmEpsilonElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->Element::GetGeometry().Area() <= 0)
        << "On EvmEpsilonElement -> " << this->Id()
        << "; Area cannot be less than or equal to 0" << std::endl;

    return 0;

    KRATOS_CATCH("");
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

template <unsigned int TDim, unsigned int TNumNodes>
std::string EvmEpsilonElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EvmEpsilonElement #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EvmEpsilonElement #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
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

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::load(Serializer& rSerializer)
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

template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                EvmEpsilonElement<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const EvmEpsilonElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

template class EvmEpsilonElement<2, 3>;
template class EvmEpsilonElement<3, 4>;
template class EvmEpsilonElement<2, 4>;
template class EvmEpsilonElement<3, 8>;

} // namespace Kratos.
