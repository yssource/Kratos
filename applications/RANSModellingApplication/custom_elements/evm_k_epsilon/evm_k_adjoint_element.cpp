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
#include "evm_k_adjoint_element.h"

#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
#include "rans_modelling_application_variables.h"

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
EvmKAdjointElement<TDim, TNumNodes>::EvmKAdjointElement(IndexType NewId)
    : StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData>(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmKAdjointElement<TDim, TNumNodes>::EvmKAdjointElement(IndexType NewId,
                                                        const NodesArrayType& ThisNodes)
    : StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData>(
          NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmKAdjointElement<TDim, TNumNodes>::EvmKAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry)
    : StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData>(
          NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmKAdjointElement<TDim, TNumNodes>::EvmKAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry,
                                                        PropertiesType::Pointer pProperties)
    : StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData>(
          NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmKAdjointElement<TDim, TNumNodes>::EvmKAdjointElement(EvmKAdjointElement<TDim, TNumNodes> const& rOther)
    : StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData>(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
EvmKAdjointElement<TDim, TNumNodes>::~EvmKAdjointElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template <unsigned int TDim, unsigned int TNumNodes>
EvmKAdjointElement<TDim, TNumNodes>& EvmKAdjointElement<TDim, TNumNodes>::operator=(
    EvmKAdjointElement<TDim, TNumNodes> const& rOther)
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
Element::Pointer EvmKAdjointElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                             NodesArrayType const& ThisNodes,
                                                             PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmKAdjointElement>(
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
Element::Pointer EvmKAdjointElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                             GeometryType::Pointer pGeom,
                                                             PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmKAdjointElement>(NewId, pGeom, pProperties);
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
Element::Pointer EvmKAdjointElement<TDim, TNumNodes>::Clone(IndexType NewId,
                                                            NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_shared<EvmKAdjointElement>(
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
void EvmKAdjointElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                           ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rResult[i] = Element::GetGeometry()[i].GetDof(RANS_ADJOINT_SCALAR_1).EquationId();
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TNumNodes)
        rElementalDofList.resize(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(RANS_ADJOINT_SCALAR_1);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    // TODO:
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::GetFirstDerivativesVector(VectorType& rValues, int Step)
{
    // TODO:
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::GetSecondDerivativesVector(VectorType& rValues, int Step)
{
    // TODO:
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod EvmKAdjointElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
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
int EvmKAdjointElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    KRATOS_CHECK_VARIABLE_KEY(TURBULENCE_RANS_C_MU);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY_SIGMA);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
    KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
    KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
    KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(RANS_AUXILIARY_VARIABLE_1);
    KRATOS_CHECK_VARIABLE_KEY(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE);
    KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS_VELOCITY_DERIVATIVES);
    KRATOS_CHECK_VARIABLE_KEY(RANS_ADJOINT_SCALAR_1);

    KRATOS_ERROR_IF(this->Id() < 1) << "EvmKAdjointElement"
                                       "found with Id 0 "
                                       "or negative"
                                    << std::endl;

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0)
        << "On EvmKAdjointElement -> " << this->Id()
        << "; Area cannot be less than or equal to 0" << std::endl;

    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
        NodeType& r_node = this->GetGeometry()[iNode];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_ADJOINT_SCALAR_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_ADJOINT_SCALAR_1, r_node);
    }

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
std::string EvmKAdjointElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EvmKAdjointElement #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EvmKAdjointElement #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
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
const Variable<double>& EvmKAdjointElement<TDim, TNumNodes>::GetPrimalVariable() const
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EvmKAdjointElement<TDim, TNumNodes>::GetPrimalVariableRelaxedRate() const
{
    return RANS_AUXILIARY_VARIABLE_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateConvectionDiffusionReactionAdjointData(
    EvmKAdjointElementData& rData,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rData.ShapeFunctionDerivatives = rShapeFunctionDerivatives;
    rData.ShapeFunctions = rShapeFunctions;

    const double& c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double& tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];

    const double& nu_t = this->EvaluateInPoint(TURBULENT_VISCOSITY, rShapeFunctions);
    const double& tke = this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, rShapeFunctions);
    const double& nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, rShapeFunctions);
    const double& wall_distance = this->EvaluateInPoint(DISTANCE, rShapeFunctions);
    const double& y_plus = this->EvaluateInPoint(RANS_Y_PLUS, rShapeFunctions);
    const double& f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
    const double& gamma = EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);

    rData.TurbulentKinematicViscosity = nu_t;
    rData.TurbulentKineticEnergy = tke;
    rData.KinematicViscosity = nu;
    rData.WallDistance = wall_distance;
    rData.Gamma = gamma;
    rData.EffectiveKinematicViscosity = nu + nu_t / tke_sigma;
    rData.Fmu = f_mu;
    rData.y_plus = y_plus;

    RansVariableUtils rans_variable_utils;

    rans_variable_utils.GetNodalArray(rData.NodalTurbulentKineticEnergy, *this,
                                      TURBULENT_KINETIC_ENERGY);
    rans_variable_utils.GetNodalArray(rData.NodalTurbulentEnergyDissipationRate,
                                      *this, TURBULENT_ENERGY_DISSIPATION_RATE);
    rans_variable_utils.GetNodalArray(rData.NodalYPlus, *this, RANS_Y_PLUS);

    std::size_t number_of_nodes = rData.NodalYPlus.size();

    if (rData.NodalFmu.size() != number_of_nodes)
        rData.NodalFmu.resize(rData.NodalYPlus.size());

    if (rData.NodalVelocity.size1() != TNumNodes || rData.NodalVelocity.size2() != TDim)
        rData.NodalVelocity.resize(TNumNodes, TDim);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rData.NodalFmu[i_node] =
            EvmKepsilonModelUtilities::CalculateFmu(rData.NodalYPlus[i_node]);
        const array_1d<double, 3>& rVelocity =
            this->GetGeometry()[i_node].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
        {
            rData.NodalVelocity(i_node, i_dim) = rVelocity[i_dim];
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateRelaxedScalarRate(
    const EvmKAdjointElementData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
{
    return this->EvaluateInPoint(RANS_AUXILIARY_VARIABLE_1, rCurrentData.ShapeFunctions);
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateEffectiveKinematicViscosity(
    const EvmKAdjointElementData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
{
    return rCurrentData.EffectiveKinematicViscosity;
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateReactionTerm(
    const EvmKAdjointElementData& rData, const ProcessInfo& rCurrentProcessInfo) const
{
    return EvmKepsilonModelUtilities::CalculateReactionTurbulentKineticEnergy(
        rData.KinematicViscosity, rData.WallDistance, rData.Gamma);
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateSourceTerm(
    const EvmKAdjointElementData& rData, const ProcessInfo& rCurrentProcessInfo) const
{
    BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
    this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rData.ShapeFunctionDerivatives);

    const double tke_production = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        velocity_gradient_matrix, rData.TurbulentKinematicViscosity, rData.TurbulentKineticEnergy);

    return tke_production;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityScalarDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const EvmKAdjointElementData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);

        noalias(rOutput) = rOutput / tke_sigma;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);

        noalias(rOutput) = rOutput / tke_sigma;
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in EvmKAdjointElement::CalculateEffectiveKinematicViscosityDerivatives method.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateReactionTermDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const EvmKAdjointElementData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);

        Vector theta_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateThetaTKESensitivity(
            theta_sensitivities, c_mu, rCurrentData.Fmu, rCurrentData.TurbulentKineticEnergy,
            rCurrentData.TurbulentKinematicViscosity, rOutput, rCurrentData.ShapeFunctions);

        noalias(rOutput) = theta_sensitivities;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);

        Vector theta_sensitivities(rOutput.size());
        EvmKepsilonModelAdjointUtilities::CalculateThetaEpsilonSensitivity(
            theta_sensitivities, c_mu, rCurrentData.Fmu, rCurrentData.TurbulentKineticEnergy,
            rCurrentData.TurbulentKinematicViscosity, rOutput);

        noalias(rOutput) = theta_sensitivities;
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in EvmKAdjointElement::CalculateReactionTermDerivatives method.";
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateSourceTermDerivatives(
    Vector& rOutput,
    const Variable<double>& rDerivativeVariable,
    const EvmKAdjointElementData& rCurrentData,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);

        BoundedMatrix<double, TDim, TDim> velocity_gradient;
        this->CalculateGradient(velocity_gradient, VELOCITY,
                                rCurrentData.ShapeFunctionDerivatives);

        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<TDim>(
            rOutput, rOutput, velocity_gradient);
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
            rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
            rCurrentData.NodalTurbulentEnergyDissipationRate, rCurrentData.NodalFmu);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            rOutput, rOutput, rCurrentData.ShapeFunctions);

        BoundedMatrix<double, TDim, TDim> velocity_gradient;
        this->CalculateGradient(velocity_gradient, VELOCITY,
                                rCurrentData.ShapeFunctionDerivatives);

        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<TDim>(
            rOutput, rOutput, velocity_gradient);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name() << " used in EvmKAdjointElement::CalculateSourceTermDerivatives method.";
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::Calculate(const Variable<Matrix>& rVariable,
                                                    Matrix& Output,
                                                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    {
        this->CalculateElementTotalResidualScalarDerivatives(
            Output, TURBULENT_ENERGY_DISSIPATION_RATE, rCurrentProcessInfo);
    }
    else
    {
        BaseType::Calculate(rVariable, Output, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityVelocityDerivatives(
    Matrix& rOutput, const EvmKAdjointElementData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
{
    const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    const Matrix& r_y_plus_sensitivities_values =
        this->GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    Matrix f_mu_sensitivities;
    EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
        f_mu_sensitivities, rCurrentData.NodalYPlus, r_y_plus_sensitivities_values);
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
        rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
        rCurrentData.NodalTurbulentEnergyDissipationRate, f_mu_sensitivities);
    EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
        rOutput, rOutput, rCurrentData.ShapeFunctions);

    noalias(rOutput) = rOutput / tke_sigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateReactionTermVelocityDerivatives(
    Matrix& rOutput, const EvmKAdjointElementData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
{
    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    const Matrix& r_y_plus_sensitivities_values =
        this->GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    Matrix f_mu_sensitivities;

    EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
        f_mu_sensitivities, rCurrentData.NodalYPlus, r_y_plus_sensitivities_values);
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
        rOutput, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
        rCurrentData.NodalTurbulentEnergyDissipationRate, f_mu_sensitivities);
    EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
        rOutput, rOutput, rCurrentData.ShapeFunctions);

    EvmKepsilonModelAdjointUtilities::CalculateGaussFmuVectorSensitivities(
        f_mu_sensitivities, rCurrentData.y_plus, r_y_plus_sensitivities_values,
        rCurrentData.ShapeFunctions);

    Matrix theta_sensitivities(rOutput.size1(), rOutput.size2());
    EvmKepsilonModelAdjointUtilities::CalculateThetaVelocitySensitivity(
        theta_sensitivities, c_mu, rCurrentData.Fmu, rCurrentData.TurbulentKineticEnergy,
        rCurrentData.TurbulentKinematicViscosity, f_mu_sensitivities, rOutput);

    noalias(rOutput) = theta_sensitivities;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::CalculateSourceTermVelocityDerivatives(
    Matrix& rOutput, const EvmKAdjointElementData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
{
    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    const Matrix& r_y_plus_sensitivities_values =
        this->GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    Matrix f_mu_sensitivities;

    EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
        f_mu_sensitivities, rCurrentData.NodalYPlus, r_y_plus_sensitivities_values);

    Matrix nu_t_sensitivities;
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
        nu_t_sensitivities, c_mu, rCurrentData.NodalTurbulentKineticEnergy,
        rCurrentData.NodalTurbulentEnergyDissipationRate, f_mu_sensitivities);
    EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
        nu_t_sensitivities, nu_t_sensitivities, rCurrentData.ShapeFunctions);

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    this->CalculateGradient(velocity_gradient, VELOCITY, rCurrentData.ShapeFunctionDerivatives);

    EvmKepsilonModelAdjointUtilities::CalculateProductionVelocitySensitivities<TDim>(
        rOutput, rCurrentData.TurbulentKinematicViscosity, nu_t_sensitivities,
        velocity_gradient, rCurrentData.ShapeFunctionDerivatives);
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityShapeSensitivity(
    const EvmKAdjointElementData& rCurrentData,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateReactionTermShapeSensitivity(
    const EvmKAdjointElementData& rCurrentData,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double EvmKAdjointElement<TDim, TNumNodes>::CalculateSourceTermShapeSensitivity(
    const EvmKAdjointElementData& rCurrentData,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    double value;

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    this->CalculateGradient(velocity_gradient, VELOCITY, rCurrentData.ShapeFunctionDerivatives);

    EvmKepsilonModelAdjointUtilities::CalculateProductionShapeSensitivities<TDim>(
        value, rCurrentData.TurbulentKinematicViscosity, 0.0, rCurrentData.NodalVelocity,
        velocity_gradient, rCurrentData.ShapeFunctionDerivatives, rDN_Dx_deriv);

    return value;
}

///@}
///@name Serialization
///@{

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes>
void EvmKAdjointElement<TDim, TNumNodes>::load(Serializer& rSerializer)
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
                                EvmKAdjointElement<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const EvmKAdjointElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

template class EvmKAdjointElement<2, 3>;
template class EvmKAdjointElement<3, 4>;
template class EvmKAdjointElement<2, 4>;
template class EvmKAdjointElement<3, 8>;

} // namespace Kratos.
