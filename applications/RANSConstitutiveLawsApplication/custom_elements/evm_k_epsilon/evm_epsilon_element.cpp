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
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::GetFirstDerivativesVector(VectorType& rValues, int Step)
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

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::GetSecondDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
            TURBULENT_ENERGY_DISSIPATION_RATE_2, Step);
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
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
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
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
        this->GetGeometryParameterDerivatives();
    const unsigned int num_gauss_points = gauss_weights.size();

    const double c1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    const double c2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];
    const int rans_time_step = rCurrentProcessInfo[RANS_TIME_STEP];
    const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
        Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                           r_parameter_derivatives_g.size2());
        noalias(contravariant_metric_tensor) =
            prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

        const array_1d<double, 3> velocity =
            this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t = this->EvaluateInPoint(
            TURBULENT_VISCOSITY, gauss_shape_functions, rans_time_step);
        const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;

        const double tke =
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon = this->EvaluateInPoint(
            TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
        const double gamma =
            EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
        const double f2 = EvmKepsilonModelUtilities::CalculateF2(nu_t, nu, epsilon);
        const double wall_distance = this->EvaluateInPoint(DISTANCE, gauss_shape_functions);
        const double reaction =
            this->CalculateReactionTerm(nu, y_plus, wall_distance, f2, c2, gamma);

        double tau, element_length;
        EvmKepsilonModelUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor, reaction,
            effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

        const double production =
            this->CalculateSourceTerm(nu_t, tke, c1, gamma, r_shape_derivatives);

        const double s = std::abs(reaction);

        for (unsigned int a = 0; a < TNumNodes; ++a)
        {
            double value = 0.0;

            value += gauss_shape_functions[a] * production;

            // Add SUPG stabilization terms
            value += (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                     tau * production;

            rRightHandSideVector[a] += gauss_weights[g] * value;
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
    KRATOS_TRY

    BaseType::CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
        this->GetGeometryParameterDerivatives();
    const unsigned int num_gauss_points = gauss_weights.size();

    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double c2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];
    const int rans_time_step = rCurrentProcessInfo[RANS_TIME_STEP];
    const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
        Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                           r_parameter_derivatives_g.size2());
        noalias(contravariant_metric_tensor) =
            prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

        const double mass = gauss_weights[g] / TNumNodes;
        this->AddLumpedMassMatrix(rMassMatrix, mass);

        const array_1d<double, 3>& velocity =
            this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t = this->EvaluateInPoint(
            TURBULENT_VISCOSITY, gauss_shape_functions, rans_time_step);
        const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;
        const double tke =
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon = this->EvaluateInPoint(
            TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double f2 = EvmKepsilonModelUtilities::CalculateF2(nu_t, nu, epsilon);
        const double wall_distance = this->EvaluateInPoint(DISTANCE, gauss_shape_functions);
        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
        const double gamma =
            EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
        const double reaction =
            this->CalculateReactionTerm(nu, y_plus, wall_distance, f2, c2, gamma);

        double tau, element_length;
        EvmKepsilonModelUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor, reaction,
            effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double s = std::abs(reaction);

        // Add mass stabilization terms
        for (unsigned int i = 0; i < TNumNodes; ++i)
            for (unsigned int j = 0; j < TNumNodes; ++j)
                rMassMatrix(i, j) +=
                    gauss_weights[g] * tau *
                    (velocity_convective_terms[i] + s * gauss_shape_functions[i]) *
                    gauss_shape_functions[j];
    }

    KRATOS_CATCH("");
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
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
        this->GetGeometryParameterDerivatives();
    const unsigned int num_gauss_points = gauss_weights.size();

    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    const double c1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    const double c2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];
    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    const double positivity_factor = rCurrentProcessInfo[RANS_STABILIZATION_MULTIPLIER];
    const int rans_time_step = rCurrentProcessInfo[RANS_TIME_STEP];
    const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma = rCurrentProcessInfo[NEWMARK_GAMMA];

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
        Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                           r_parameter_derivatives_g.size2());
        noalias(contravariant_metric_tensor) =
            prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

        const double epsilon = this->EvaluateInPoint(
            TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double tke =
            this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double nu_t = this->EvaluateInPoint(
            TURBULENT_VISCOSITY, gauss_shape_functions, rans_time_step);
        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;

        const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
        const double gamma =
            EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);

        array_1d<double, 3> velocity =
            this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double velocity_magnitude = norm_2(velocity);
        const double f2 = EvmKepsilonModelUtilities::CalculateF2(tke, nu, epsilon);
        const double wall_distance = this->EvaluateInPoint(DISTANCE, gauss_shape_functions);
        const double reaction =
            this->CalculateReactionTerm(nu, y_plus, wall_distance, f2, c2, gamma);

        double tau, element_length;
        EvmKepsilonModelUtilities::CalculateStabilizationTau(
            tau, element_length, velocity, contravariant_metric_tensor, reaction,
            effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time);

        // Calculate residual for cross wind dissipation coefficient
        array_1d<double, 3> epsilon_gradient;
        this->CalculateGradient(epsilon_gradient, TURBULENT_ENERGY_DISSIPATION_RATE,
                                r_shape_derivatives);
        const double epsilon_gradient_norm = norm_2(epsilon_gradient);

        double cross_wind_diffusion{0.0}, stream_line_diffusion{0.0};
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

        Vector nodal_epsilon;
        this->GetValuesVector(nodal_epsilon);

        if (epsilon_gradient_norm > std::numeric_limits<double>::epsilon() &&
            velocity_magnitude_square > std::numeric_limits<double>::epsilon())
        {
            const double relaxed_epsilon_acceleration = this->EvaluateInPoint(
                RANS_AUXILIARY_VARIABLE_2, gauss_shape_functions);
            const double source =
                this->CalculateSourceTerm(nu_t, tke, c1, gamma, r_shape_derivatives);

            double residual = relaxed_epsilon_acceleration;
            residual += inner_prod(velocity_convective_terms, nodal_epsilon);
            residual += reaction * inner_prod(gauss_shape_functions, nodal_epsilon);
            residual -= source;
            residual = std::abs(residual);
            residual /= epsilon_gradient_norm;

            double chi, k1, k2;
            EvmKepsilonModelUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                reaction, bossak_alpha, bossak_gamma, delta_time, element_length);

            stream_line_diffusion = residual * chi * k1 / velocity_magnitude_square;
            cross_wind_diffusion = residual * chi * k2 / velocity_magnitude_square;
        }

        const double s = std::abs(reaction);

        stream_line_diffusion *= positivity_factor;
        cross_wind_diffusion *= positivity_factor;

        for (unsigned int a = 0; a < TNumNodes; a++)
        {
            for (unsigned int b = 0; b < TNumNodes; b++)
            {
                double dNa_dNb = 0.0;
                for (unsigned int i = 0; i < TDim; i++)
                    dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);

                double value = 0.0;

                value += gauss_shape_functions[a] * velocity_convective_terms[b];
                value += gauss_shape_functions[a] * reaction *
                         gauss_shape_functions[b]; // * positive_values_list[b];
                value += effective_kinematic_viscosity * dNa_dNb;

                // Adding SUPG stabilization terms
                value += tau *
                         (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                         velocity_convective_terms[b];
                value +=
                    tau * (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                    reaction * gauss_shape_functions[b]; // * positive_values_list[b];

                // Adding cross wind dissipation
                value += cross_wind_diffusion * dNa_dNb * velocity_magnitude_square;
                value -= cross_wind_diffusion * velocity_convective_terms[a] *
                         velocity_convective_terms[b];

                // Adding stream line dissipation
                value += stream_line_diffusion * velocity_convective_terms[a] *
                         velocity_convective_terms[b];

                rDampingMatrix(a, b) += gauss_weights[g] * value;
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);

    // Now calculate an additional contribution to the residual: r -= rDampingMatrix * (u,p)
    VectorType U = ZeroVector(TNumNodes);
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        U[iNode] = this->GetGeometry()[iNode].FastGetSolutionStepValue(
            TURBULENT_ENERGY_DISSIPATION_RATE);

    noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmEpsilonElement<TDim, TNumNodes>::Calculate(const Variable<double>& rVariable,
                                                   double& Output,
                                                   const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == RESIDUAL)
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        Output = 0.0;

        // const double epsilon_sigma =
        //     rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
        const double c1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
        const double c2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
        // const double delta_time = rCurrentProcessInfo[DELTA_TIME];
        // const double positivity_factor = rCurrentProcessInfo[RANS_STABILIZATION_MULTIPLIER];
        const int rans_time_step = rCurrentProcessInfo[RANS_TIME_STEP];

        Vector nodal_epsilon;
        this->GetValuesVector(nodal_epsilon);

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const double epsilon = this->EvaluateInPoint(
                TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
            const double tke =
                this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
            const double nu_t = this->EvaluateInPoint(
                TURBULENT_VISCOSITY, gauss_shape_functions, rans_time_step);
            const double nu =
                this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
            // const double effective_kinematic_viscosity = nu + nu_t / epsilon_sigma;

            const double y_plus = this->EvaluateInPoint(RANS_Y_PLUS, gauss_shape_functions);
            const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
            const double gamma =
                EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);

            array_1d<double, 3> velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            // const double velocity_magnitude = norm_2(velocity);
            const double f2 = EvmKepsilonModelUtilities::CalculateF2(tke, nu, epsilon);
            const double wall_distance =
                this->EvaluateInPoint(DISTANCE, gauss_shape_functions);
            const double reaction =
                this->CalculateReactionTerm(nu, y_plus, wall_distance, f2, c2, gamma);

            const double relaxed_epsilon_acceleration = this->EvaluateInPoint(
                RANS_AUXILIARY_VARIABLE_2, gauss_shape_functions);
            const double source =
                this->CalculateSourceTerm(nu_t, tke, c1, gamma, r_shape_derivatives);

            Output += relaxed_epsilon_acceleration;
            Output += inner_prod(velocity_convective_terms, nodal_epsilon);
            Output += reaction * inner_prod(gauss_shape_functions, nodal_epsilon);
            Output -= source;
            Output = std::abs(Output);
        }
    }
    else
    {
        KRATOS_ERROR
            << "Unsupported variable used in EvmKElement Calculate method"
            << rVariable << "\n.";
    }
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

template <unsigned int TDim, unsigned int TNumNodes>
double EvmEpsilonElement<TDim, TNumNodes>::CalculateReactionTerm(const double kinematic_viscosity,
                                                                 const double y_plus,
                                                                 const double wall_distance,
                                                                 const double f2,
                                                                 const double c2,
                                                                 const double gamma) const
{
    return c2 * f2 * gamma + 2.0 * kinematic_viscosity * std::exp(-0.5 * y_plus) /
                                 std::pow(wall_distance, 2);
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmEpsilonElement<TDim, TNumNodes>::CalculateSourceTerm(const double turbulent_kinematic_viscosity,
                                                               const double turbulent_kinetic_energy,
                                                               const double c1,
                                                               const double gamma,
                                                               const Matrix& rShapeDerivatives) const
{
    BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
    this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rShapeDerivatives);
    double production = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        velocity_gradient_matrix, turbulent_kinematic_viscosity, turbulent_kinetic_energy);

    production *= (c1 * gamma);

    return production;
}

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
