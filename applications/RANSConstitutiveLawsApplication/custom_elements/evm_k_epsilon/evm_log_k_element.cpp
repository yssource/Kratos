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
#include "evm_log_k_element.h"

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
EvmLogKElement<TDim, TNumNodes>::EvmLogKElement(IndexType NewId)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmLogKElement<TDim, TNumNodes>::EvmLogKElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmLogKElement<TDim, TNumNodes>::EvmLogKElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmLogKElement<TDim, TNumNodes>::EvmLogKElement(IndexType NewId,
                                                GeometryType::Pointer pGeometry,
                                                PropertiesType::Pointer pProperties)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmLogKElement<TDim, TNumNodes>::EvmLogKElement(EvmLogKElement<TDim, TNumNodes> const& rOther)
    : RANSConstitutiveElement<TDim, TNumNodes, 1>(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmLogKElement<TDim, TNumNodes>::~EvmLogKElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template <unsigned int TDim, unsigned int TNumNodes>
EvmLogKElement<TDim, TNumNodes>& EvmLogKElement<TDim, TNumNodes>::operator=(
    EvmLogKElement<TDim, TNumNodes> const& rOther)
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
Element::Pointer EvmLogKElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                         NodesArrayType const& ThisNodes,
                                                         PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmLogKElement>(
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
Element::Pointer EvmLogKElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                         GeometryType::Pointer pGeom,
                                                         PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmLogKElement>(NewId, pGeom, pProperties);
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
Element::Pointer EvmLogKElement<TDim, TNumNodes>::Clone(IndexType NewId,
                                                        NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmLogKElement>(
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
void EvmLogKElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                       ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != TLocalSize)
        rResult.resize(TLocalSize, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rResult[i] =
            Element::GetGeometry()[i].GetDof(TURBULENT_LOG_KINETIC_ENERGY).EquationId();
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TLocalSize)
        rElementalDofList.resize(TLocalSize);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rElementalDofList[i] =
            Element::GetGeometry()[i].pGetDof(TURBULENT_LOG_KINETIC_ENERGY);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::GetFirstDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] =
            rGeom[iNode].FastGetSolutionStepValue(TURBULENT_LOG_KINETIC_ENERGY, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::GetSecondDerivativesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
            TURBULENT_LOG_KINETIC_ENERGY_RATE, Step);
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
void EvmLogKElement<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                           VectorType& rRightHandSideVector,
                                                           ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
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
void EvmLogKElement<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
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
void EvmLogKElement<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
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

    const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
        Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                           r_parameter_derivatives_g.size2());
        noalias(contravariant_metric_tensor) =
            prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

        const array_1d<double, 3> effective_convective_velocity =
            this->GetConvectiveVelocity(tke_sigma, gauss_shape_functions, r_shape_derivatives);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms,
                                    effective_convective_velocity, r_shape_derivatives);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t =
            this->EvaluateInPoint(TURBULENT_VISCOSITY, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / tke_sigma;

        double tau, element_length;
        EvmLogKepsilonModelUtilities::CalculateStabilizationTau(
            tau, element_length, effective_convective_velocity,
            contravariant_metric_tensor, effective_kinematic_viscosity, delta_time);

        const double log_tke = this->EvaluateInPoint(TURBULENT_LOG_KINETIC_ENERGY, gauss_shape_functions);
        const double log_epsilon = this->EvaluateInPoint(TURBULENT_LOG_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double wall_distance = this->EvaluateInPoint(DISTANCE, gauss_shape_functions);

        const double gamma = EvmLogKepsilonModelUtilities::CalculateGamma(log_tke, log_epsilon);
        const double source = this->CalculateSourceTerm(gamma, nu, nu_t, wall_distance, log_tke, r_shape_derivatives);

        for (unsigned int a = 0; a < TNumNodes; ++a)
        {
            double value = 0.0;

            value += gauss_shape_functions[a] * source;

            // Add supg stabilization terms
            value += velocity_convective_terms[a] * tau * source;

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
void EvmLogKElement<TDim, TNumNodes>::CalculateFirstDerivativesContributions(
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
void EvmLogKElement<TDim, TNumNodes>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
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
void EvmLogKElement<TDim, TNumNodes>::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
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
void EvmLogKElement<TDim, TNumNodes>::CalculateSecondDerivativesContributions(
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
void EvmLogKElement<TDim, TNumNodes>::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                                    ProcessInfo& rCurrentProcessInfo)
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
void EvmLogKElement<TDim, TNumNodes>::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                                                    ProcessInfo& rCurrentProcessInfo)
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
void EvmLogKElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
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
    const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];

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

        const array_1d<double, 3> effective_convective_velocity =
            this->GetConvectiveVelocity(tke_sigma, gauss_shape_functions, r_shape_derivatives);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms,
                                    effective_convective_velocity, r_shape_derivatives);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t =
            this->EvaluateInPoint(TURBULENT_VISCOSITY, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / tke_sigma;

        double tau, element_length;
        EvmLogKepsilonModelUtilities::CalculateStabilizationTau(
            tau, element_length, effective_convective_velocity,
            contravariant_metric_tensor, effective_kinematic_viscosity, delta_time);

        // Add mass stabilization terms
        for (unsigned int i = 0; i < TNumNodes; ++i)
            for (unsigned int j = 0; j < TNumNodes; ++j)
                rMassMatrix(i, j) += gauss_weights[g] * tau *
                                     velocity_convective_terms[i] *
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
void EvmLogKElement<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix,
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

    const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    for (unsigned int g = 0; g < num_gauss_points; g++)
    {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

        Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                           r_parameter_derivatives_g.size2());
        noalias(contravariant_metric_tensor) =
            prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

        array_1d<double, 3> effective_convective_velocity = this->GetConvectiveVelocity(
            tke_sigma, gauss_shape_functions, r_shape_derivatives);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms,
                                    effective_convective_velocity, r_shape_derivatives);
        const double effective_velocity_magnitude_square =
            std::pow(norm_2(effective_convective_velocity), 2);

        const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, gauss_shape_functions);
        const double nu_t =
            this->EvaluateInPoint(TURBULENT_VISCOSITY, gauss_shape_functions);
        const double effective_kinematic_viscosity = nu + nu_t / tke_sigma;

        double tau, element_length;
        EvmLogKepsilonModelUtilities::CalculateStabilizationTau(
            tau, element_length, effective_convective_velocity,
            contravariant_metric_tensor, effective_kinematic_viscosity, delta_time);

        // Calculate residual for cross wind dissipation coefficient
        double cross_wind_diffusion{0.0}, stream_line_diffusion{0.0};

        array_1d<double, 3> log_tke_gradient;
        this->CalculateGradient(log_tke_gradient, TURBULENT_LOG_KINETIC_ENERGY,
                                r_shape_derivatives);
        const double log_tke_gradient_norm = norm_2(log_tke_gradient);

        if (log_tke_gradient_norm > std::numeric_limits<double>::epsilon() &&
            effective_velocity_magnitude_square > std::numeric_limits<double>::epsilon())
        {
            const double relaxed_log_tke_acceleration = this->EvaluateInPoint(
                RANS_AUXILIARY_VARIABLE_1, gauss_shape_functions);
            const double log_k = this->EvaluateInPoint(
                TURBULENT_LOG_KINETIC_ENERGY, gauss_shape_functions);
            const double log_epsilon = this->EvaluateInPoint(
                TURBULENT_LOG_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
            const double wall_distance =
                this->EvaluateInPoint(DISTANCE, gauss_shape_functions);

            const double gamma =
                EvmLogKepsilonModelUtilities::CalculateGamma(log_k, log_epsilon);
            const double source = this->CalculateSourceTerm(
                gamma, nu, nu_t, wall_distance, log_k, r_shape_derivatives);

            Vector nodal_log_tke;
            this->GetValuesVector(nodal_log_tke);

            double residual = relaxed_log_tke_acceleration;
            residual += inner_prod(velocity_convective_terms, nodal_log_tke);
            residual -= source;
            residual = std::abs(residual);
            residual /= log_tke_gradient_norm;

            double chi, k1, k2;
            EvmLogKepsilonModelUtilities::CalculateCrossWindDiffusionParameters(
                chi, k1, k2, std::sqrt(effective_velocity_magnitude_square),
                tau, effective_kinematic_viscosity, 0.0, element_length);

            stream_line_diffusion = residual * chi * k1 / effective_velocity_magnitude_square;
            cross_wind_diffusion = residual * chi * k2 / effective_velocity_magnitude_square;
        }

        for (unsigned int a = 0; a < TNumNodes; a++)
        {
            for (unsigned int b = 0; b < TNumNodes; b++)
            {
                double dNa_dNb = 0.0;
                for (unsigned int i = 0; i < TDim; i++)
                    dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);

                double value = 0.0;

                value += gauss_shape_functions[a] * velocity_convective_terms[b];
                value += effective_kinematic_viscosity * dNa_dNb;

                // Adding SUPG stabilization terms
                value += tau * (velocity_convective_terms[a]) *
                         velocity_convective_terms[b];

                // Adding cross wind dissipation
                value += cross_wind_diffusion * dNa_dNb * effective_velocity_magnitude_square;
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
void EvmLogKElement<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);

    // Now calculate an additional contribution to the residual: r -= rDampingMatrix * (u,p)
    VectorType U = ZeroVector(TNumNodes);
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        U[iNode] = this->GetGeometry()[iNode].FastGetSolutionStepValue(TURBULENT_LOG_KINETIC_ENERGY);

    noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
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
int EvmLogKElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)
        << "EvmLogKElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->Element::GetGeometry().Area() <= 0)
        << "On EvmLogKElement -> " << this->Id()
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
std::string EvmLogKElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EvmLogKElement #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EvmLogKElement #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
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
array_1d<double, 3> EvmLogKElement<TDim, TNumNodes>::GetConvectiveVelocity(
    const double turbulent_kinetic_energy_sigma,
    const Vector& rGaussShapeFunctions,
    const Matrix& rGaussShapeDerivatives) const
{
    const array_1d<double, 3> velocity =
        this->EvaluateInPoint(VELOCITY, rGaussShapeFunctions);
    array_1d<double, 3> log_k_gradient;
    this->CalculateGradient(log_k_gradient, TURBULENT_LOG_KINETIC_ENERGY, rGaussShapeDerivatives);
    const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, rGaussShapeFunctions);
    const double nu_t = this->EvaluateInPoint(TURBULENT_VISCOSITY, rGaussShapeFunctions);
    const double effective_kinematic_viscosity = nu + nu_t / turbulent_kinetic_energy_sigma;

    return velocity - effective_kinematic_viscosity * log_k_gradient;
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmLogKElement<TDim, TNumNodes>::CalculateSourceTerm(
    const double gamma,
    const double kinematic_viscosity,
    const double turbulent_kinematic_viscosity,
    const double wall_distance,
    const double log_turbulent_kinetic_energy,
    const Matrix& rShapeDerivatives) const
{
    const double reaction = 2.0 * kinematic_viscosity / std::pow(wall_distance, 2) + gamma;

    BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
    this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rShapeDerivatives);

    const double P = EvmLogKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        velocity_gradient_matrix, turbulent_kinematic_viscosity);

    return P * std::exp(-1.0 * log_turbulent_kinetic_energy) - reaction;
}

///@}
///@name Serialization
///@{

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmLogKElement<TDim, TNumNodes>::load(Serializer& rSerializer)
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
                                EvmLogKElement<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const EvmLogKElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

template class EvmLogKElement<2, 3>;
template class EvmLogKElement<3, 4>;
template class EvmLogKElement<2, 4>;
template class EvmLogKElement<3, 8>;

} // namespace Kratos.
