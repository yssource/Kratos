//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <random>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_process.h"
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_sensitivities_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/model_part.h"
#include "rans_modelling_application_variables.h"
#include "test_utilities.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
namespace Testing
{
typedef ModelPart::NodeType NodeType;

typedef ModelPart::ElementType ElementType;

typedef Geometry<NodeType> GeometryType;

typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

/**
 * Auxiliar function to generate a triangular element to be tested.
 */
namespace RansEvmKEpsilonModel
{
void AddVariablesToModelPart(ModelPart& rModelPart)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(RELAXED_ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    rModelPart.AddNodalSolutionStepVariable(RANS_ADJOINT_SCALAR_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_ADJOINT_SCALAR_2);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);
}

void InitializeProcessInfo(ModelPart& rModelPart)
{
    // Process info creation
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.1);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 0.01);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 0.09);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C1, 1.44);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C2, 1.92);
    rModelPart.GetProcessInfo().SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 1.03);
    rModelPart.GetProcessInfo().SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.3);
    rModelPart.GetProcessInfo().SetValue(BOSSAK_ALPHA, -0.03);
    rModelPart.GetProcessInfo().SetValue(NEWMARK_GAMMA, 1.25);
    rModelPart.GetProcessInfo().SetValue(WALL_SMOOTHNESS_BETA, 5.2);
    rModelPart.GetProcessInfo().SetValue(WALL_VON_KARMAN, 0.41);
    rModelPart.GetProcessInfo().SetValue(OSS_SWITCH, 0);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 3.0e-02);
}

void CreateModelPartNodes(ModelPart& rModelPart)
{
    // Element creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

    VariableUtils().AddDof<Variable<double>>(TURBULENT_KINETIC_ENERGY, rModelPart);
    VariableUtils().AddDof<Variable<double>>(TURBULENT_ENERGY_DISSIPATION_RATE, rModelPart);
    VariableUtils().AddDof<Variable<double>>(RANS_ADJOINT_SCALAR_1, rModelPart);
    VariableUtils().AddDof<Variable<double>>(RANS_ADJOINT_SCALAR_2, rModelPart);

    // VariableUtils().AddDof<Variable<array_1d<double, 3>>>(ADJOINT_FLUID_VECTOR_1, rModelPart);
    // VariableUtils().AddDof<Variable<array_1d<double, 3>>>(ADJOINT_FLUID_VECTOR_2, rModelPart);
    // VariableUtils().AddDof<Variable<array_1d<double, 3>>>(ADJOINT_FLUID_VECTOR_3, rModelPart);
    // VariableUtils().AddDof<Variable<double>>(ADJOINT_FLUID_SCALAR_1, rModelPart);
}

void CreateModelPartElements(ModelPart& rModelPart, std::string ElementName)
{
    Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);

    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
    rModelPart.CreateNewElement(ElementName,
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_nodes, p_elem_prop);
}

void InitializeNodalVariables(ModelPart& rModelPart)
{
    using namespace RansModellingApplicationTestUtilities;

    InitializeVariableWithRandomValues(rModelPart, DISTANCE, 1e-2, 1.0, 2);
    InitializeVariableWithRandomValues(rModelPart, TURBULENT_KINETIC_ENERGY, 0.1, 1.0, 2);
    InitializeVariableWithRandomValues(
        rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 5.0, 10.0, 2);
    InitializeVariableWithRandomValues(
        rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 10, 20.0, 2);
    InitializeVariableWithRandomValues(
        rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 5, 10.0, 2);
    InitializeVariableWithRandomValues(rModelPart, VELOCITY, 5, 10.0, 2);

    InitializeVariableWithRandomValues(rModelPart, PRESSURE, 5, 10.0, 2);
    InitializeVariableWithRandomValues(rModelPart, ACCELERATION, 2.0, 3.0, 2);

    InitializeVariableWithValues(rModelPart, KINEMATIC_VISCOSITY, 3e-2);
    InitializeVariableWithValues(rModelPart, DENSITY, 200.0);
}

void GenerateRansEvmKEpsilonTestModelPart(ModelPart& rModelPart, std::string ElementName)
{
    AddVariablesToModelPart(rModelPart);
    InitializeProcessInfo(rModelPart);
    CreateModelPartNodes(rModelPart);
    CreateModelPartElements(rModelPart, ElementName);
    InitializeNodalVariables(rModelPart);
}

void UpdateVariablesInModelPart(ModelPart& rModelPart)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];
    const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);


        double& nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke, epsilon, f_mu);
        r_node.FastGetSolutionStepValue(VISCOSITY) = nu + nu_t;

        const double tke_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
        const double old_tke_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 1);
        const double epsilon_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        const double old_epsilon_rate =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2, 1);
        r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) =
            (1 - bossak_alpha) * tke_rate + bossak_alpha * old_tke_rate;
        r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) =
            (1 - bossak_alpha) * epsilon_rate + bossak_alpha * old_epsilon_rate;

        const array_1d<double, 3>& acceleration =
            r_node.FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double, 3>& old_acceleration =
            r_node.FastGetSolutionStepValue(ACCELERATION, 1);

        r_node.FastGetSolutionStepValue(RELAXED_ACCELERATION) =
            acceleration * (1 - bossak_alpha) + bossak_alpha * old_acceleration;
    }
}

void CalculatePrimalQuantities(std::vector<double>& rValues,
                               const ElementType& rElement,
                               const Vector& rGaussShapeFunctions,
                               const Matrix& rGaussShapeFunctionDerivatives,
                               const ProcessInfo& rCurrentProcessInfo)
{
    RansCalculationUtilities rans_calculation_utilities;

    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    const GeometryType& r_geometry = rElement.GetGeometry();

    const double y_plus = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, RANS_Y_PLUS, rGaussShapeFunctions);
    const double tke = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_KINETIC_ENERGY, rGaussShapeFunctions);
    const double epsilon = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, rGaussShapeFunctions);
    const double nu_t = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_VISCOSITY, rGaussShapeFunctions);
    const double nu = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, KINEMATIC_VISCOSITY, rGaussShapeFunctions);

    const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
    const double Re_t = std::pow(tke, 2) / (nu * epsilon);
    const double theta = EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
    const double f2 = EvmKepsilonModelUtilities::CalculateF2(tke, nu, epsilon);

    BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
    rans_calculation_utilities.CalculateGradient<2>(
        velocity_gradient_matrix, r_geometry, VELOCITY, rGaussShapeFunctionDerivatives);
    const double P_k = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
        velocity_gradient_matrix, nu_t, tke);

    rValues.clear();
    rValues.push_back(nu_t);
    rValues.push_back(P_k);
    rValues.push_back(theta);
    rValues.push_back(Re_t);
    rValues.push_back(f2);
    rValues.push_back(f_mu);
    rValues.push_back(nu);
    rValues.push_back(epsilon);
    rValues.push_back(tke);
    rValues.push_back(y_plus);
}

void ReadNodalDataFromElement(Vector& rYPlus,
                              Vector& rTKE,
                              Vector& rEpsilon,
                              Vector& rNut,
                              Vector& rFmu,
                              const Element& rElement)
{
    const auto& r_geometry = rElement.GetGeometry();
    std::size_t number_of_nodes = r_geometry.PointsNumber();

    RansVariableUtils rans_variable_utils;

    rans_variable_utils.GetNodalArray(rTKE, rElement, TURBULENT_KINETIC_ENERGY);
    rans_variable_utils.GetNodalArray(rEpsilon, rElement, TURBULENT_ENERGY_DISSIPATION_RATE);
    rans_variable_utils.GetNodalArray(rYPlus, rElement, RANS_Y_PLUS);
    rans_variable_utils.GetNodalArray(rNut, rElement, TURBULENT_VISCOSITY);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        rFmu[i_node] = EvmKepsilonModelUtilities::CalculateFmu(rYPlus[i_node]);
}
} // namespace RansEvmKEpsilonModel

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussTKESensitivities, RANSEvModelsKEpsilonGaussMatrices)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("RansSensitivities");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelProcess y_plus_model_process(r_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivities = [](std::vector<Vector>& rValues,
                                      const ElementType& rElement,
                                      const Vector& rGaussShapeFunctions,
                                      const Matrix& rGaussShapeFunctionDerivatives,
                                      const ProcessInfo& rCurrentProcessInfo) {
        RansCalculationUtilities rans_calculation_utilities;

        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        const GeometryType& r_geometry = rElement.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Vector nodal_y_plus(number_of_nodes);
        Vector nodal_tke(number_of_nodes);
        Vector nodal_epsilon(number_of_nodes);
        Vector nodal_nu_t(number_of_nodes);
        Vector nodal_f_mu(number_of_nodes);

        RansEvmKEpsilonModel::ReadNodalDataFromElement(
            nodal_y_plus, nodal_tke, nodal_epsilon, nodal_nu_t, nodal_f_mu, rElement);

        std::vector<double> primal_quantities;
        RansEvmKEpsilonModel::CalculatePrimalQuantities(
            primal_quantities, rElement, rGaussShapeFunctions,
            rGaussShapeFunctionDerivatives, rCurrentProcessInfo);

        const double nu_t = primal_quantities[0];
        const double Re_t = primal_quantities[3];
        const double f_mu = primal_quantities[5];
        const double nu = primal_quantities[6];
        const double epsilon = primal_quantities[7];
        const double tke = primal_quantities[8];

        Vector nodal_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
            nodal_nu_t_sensitivities, c_mu, nodal_tke, nodal_epsilon, nodal_f_mu);

        Vector gauss_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_nu_t_sensitivities, nodal_nu_t_sensitivities, rGaussShapeFunctions);

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix, r_geometry, VELOCITY, rGaussShapeFunctionDerivatives);

        Vector gauss_production_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<2>(
            gauss_production_sensitivities, gauss_nu_t_sensitivities, velocity_gradient_matrix);

        Vector gauss_re_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateTurbulentReynoldsNumberTKESensitivity(
            gauss_re_t_sensitivities, tke, epsilon, nu, rGaussShapeFunctions);

        Vector gauss_theta_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateThetaTKESensitivity(
            gauss_theta_sensitivities, c_mu, f_mu, tke, nu_t,
            gauss_nu_t_sensitivities, rGaussShapeFunctions);

        Vector gauss_f2_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateF2ScalarSensitivity(
            gauss_f2_sensitivities, epsilon, Re_t, gauss_re_t_sensitivities);

        rValues.clear();
        rValues.push_back(gauss_nu_t_sensitivities);
        rValues.push_back(gauss_production_sensitivities);
        rValues.push_back(gauss_theta_sensitivities);
        rValues.push_back(gauss_re_t_sensitivities);
        rValues.push_back(gauss_f2_sensitivities);
    };

    RansModellingApplicationTestUtilities::RunGaussPointScalarSensitivityTest(
        r_model_part, y_plus_model_process,
        RansEvmKEpsilonModel::CalculatePrimalQuantities, calculate_sensitivities,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, perturb_variable, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussEpsilonSensitivities, RANSEvModelsKEpsilonGaussMatrices)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("RansSensitivities");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelProcess y_plus_model_process(r_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivities = [](std::vector<Vector>& rValues,
                                      const ElementType& rElement,
                                      const Vector& rGaussShapeFunctions,
                                      const Matrix& rGaussShapeFunctionDerivatives,
                                      const ProcessInfo& rCurrentProcessInfo) {
        RansCalculationUtilities rans_calculation_utilities;

        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        const GeometryType& r_geometry = rElement.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Vector nodal_y_plus(number_of_nodes);
        Vector nodal_tke(number_of_nodes);
        Vector nodal_epsilon(number_of_nodes);
        Vector nodal_nu_t(number_of_nodes);
        Vector nodal_f_mu(number_of_nodes);

        RansEvmKEpsilonModel::ReadNodalDataFromElement(
            nodal_y_plus, nodal_tke, nodal_epsilon, nodal_nu_t, nodal_f_mu, rElement);

        std::vector<double> primal_quantities;
        RansEvmKEpsilonModel::CalculatePrimalQuantities(
            primal_quantities, rElement, rGaussShapeFunctions,
            rGaussShapeFunctionDerivatives, rCurrentProcessInfo);

        const double nu_t = primal_quantities[0];
        const double Re_t = primal_quantities[3];
        const double f_mu = primal_quantities[5];
        const double nu = primal_quantities[6];
        const double epsilon = primal_quantities[7];
        const double tke = primal_quantities[8];

        Vector nodal_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
            nodal_nu_t_sensitivities, c_mu, nodal_tke, nodal_epsilon, nodal_f_mu);

        Vector gauss_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_nu_t_sensitivities, nodal_nu_t_sensitivities, rGaussShapeFunctions);

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix, r_geometry, VELOCITY, rGaussShapeFunctionDerivatives);

        Vector gauss_production_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<2>(
            gauss_production_sensitivities, gauss_nu_t_sensitivities, velocity_gradient_matrix);

        Vector gauss_theta_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateThetaEpsilonSensitivity(
            gauss_theta_sensitivities, c_mu, f_mu, tke, nu_t, gauss_nu_t_sensitivities);

        Vector gauss_re_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateTurbulentReynoldsNumberEpsilonSensitivity(
            gauss_re_t_sensitivities, tke, epsilon, nu, rGaussShapeFunctions);

        Vector gauss_f2_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateF2ScalarSensitivity(
            gauss_f2_sensitivities, epsilon, Re_t, gauss_re_t_sensitivities);

        rValues.clear();
        rValues.push_back(gauss_nu_t_sensitivities);
        rValues.push_back(gauss_production_sensitivities);
        rValues.push_back(gauss_theta_sensitivities);
        rValues.push_back(gauss_re_t_sensitivities);
        rValues.push_back(gauss_f2_sensitivities);
    };

    RansModellingApplicationTestUtilities::RunGaussPointScalarSensitivityTest(
        r_model_part, y_plus_model_process,
        RansEvmKEpsilonModel::CalculatePrimalQuantities, calculate_sensitivities,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, perturb_variable, 1e-7, 3e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussVelocitySensitivities, RANSEvModelsKEpsilonGaussMatrices)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("RansSensitivities");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_model_sensitivities_process(
        r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess y_plus_model_process(r_model_part, empty_parameters);

    y_plus_model_process.Check();
    y_plus_model_process.Execute();

    y_plus_model_sensitivities_process.Check();
    y_plus_model_sensitivities_process.Execute();

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_vector = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_vector[Dim];
    };

    auto calculate_sensitivities = [](std::vector<Matrix>& rValues,
                                      const ElementType& rElement,
                                      const Vector& rGaussShapeFunctions,
                                      const Matrix& rGaussShapeFunctionDerivatives,
                                      const ProcessInfo& rCurrentProcessInfo) {
        RansCalculationUtilities rans_calculation_utilities;

        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
        const int domain_size = rCurrentProcessInfo[DOMAIN_SIZE];

        const GeometryType& r_geometry = rElement.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Vector nodal_y_plus(number_of_nodes);
        Vector nodal_tke(number_of_nodes);
        Vector nodal_epsilon(number_of_nodes);
        Vector nodal_nu_t(number_of_nodes);
        Vector nodal_f_mu(number_of_nodes);

        RansEvmKEpsilonModel::ReadNodalDataFromElement(
            nodal_y_plus, nodal_tke, nodal_epsilon, nodal_nu_t, nodal_f_mu, rElement);

        std::vector<double> primal_quantities;
        RansEvmKEpsilonModel::CalculatePrimalQuantities(
            primal_quantities, rElement, rGaussShapeFunctions,
            rGaussShapeFunctionDerivatives, rCurrentProcessInfo);

        const double nu_t = primal_quantities[0];
        const double f_mu = primal_quantities[5];
        const double tke = primal_quantities[8];
        const double y_plus = primal_quantities[9];

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix, r_geometry, VELOCITY, rGaussShapeFunctionDerivatives);

        const Matrix& y_plus_sensitivities =
            rElement.GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

        Matrix nodal_f_mu_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
            nodal_f_mu_sensitivities, nodal_y_plus, y_plus_sensitivities);

        Matrix nodal_nu_t_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
            nodal_nu_t_sensitivities, c_mu, nodal_tke, nodal_epsilon, nodal_f_mu_sensitivities);

        Matrix gauss_f_mu_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussFmuVectorSensitivities(
            gauss_f_mu_velocity_sensitivities, y_plus, y_plus_sensitivities, rGaussShapeFunctions);

        Matrix gauss_nu_t_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_nu_t_sensitivities, nodal_nu_t_sensitivities, rGaussShapeFunctions);

        Matrix gauss_y_plus_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_y_plus_velocity_sensitivities, y_plus_sensitivities, rGaussShapeFunctions);

        Matrix gauss_production_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateProductionVelocitySensitivities<2>(
            gauss_production_sensitivities, nu_t, gauss_nu_t_sensitivities,
            velocity_gradient_matrix, rGaussShapeFunctionDerivatives);

        Matrix gauss_theta_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateThetaVelocitySensitivity(
            gauss_theta_sensitivities, c_mu, f_mu, tke, nu_t,
            gauss_f_mu_velocity_sensitivities, gauss_nu_t_sensitivities);

        rValues.clear();
        rValues.push_back(gauss_nu_t_sensitivities);
        rValues.push_back(gauss_production_sensitivities);
        rValues.push_back(gauss_theta_sensitivities);
    };

    RansModellingApplicationTestUtilities::RunGaussPointVectorSensitivityTest(
        r_model_part, y_plus_model_process,
        RansEvmKEpsilonModel::CalculatePrimalQuantities, calculate_sensitivities,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart, perturb_variable, 1e-7, 3e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementTKEFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementEpsilonFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementTKESecondDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementVelocityDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_VELOCITY_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementShapeSensitivity, RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMK2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonElementEpsilonFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMEPSILON2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMEpsilonAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonElementTKEFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMEPSILON2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMEpsilonAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-8, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonElementEpsilonSecondDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMEPSILON2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMEpsilonAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonElementShapeSensitivity,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMEPSILON2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMEpsilonAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmEpsilonElementVelocityDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "RANSEVMEPSILON2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMEpsilonAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_VELOCITY_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonNodalTKESensitivities, RANSEvModelsKEpsilonNodalMatrices)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("RansSensitivities");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelProcess y_plus_model_process(r_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivities = [](std::vector<Vector>& rValues,
                                      const ElementType& rElement,
                                      const ProcessInfo& rCurrentProcessInfo) {
        RansCalculationUtilities rans_calculation_utilities;

        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        const GeometryType& r_geometry = rElement.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Vector nodal_y_plus(number_of_nodes);
        Vector nodal_tke(number_of_nodes);
        Vector nodal_epsilon(number_of_nodes);
        Vector nodal_nu_t(number_of_nodes);
        Vector nodal_f_mu(number_of_nodes);

        RansEvmKEpsilonModel::ReadNodalDataFromElement(
            nodal_y_plus, nodal_tke, nodal_epsilon, nodal_nu_t, nodal_f_mu, rElement);

        Vector nodal_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
            nodal_nu_t_sensitivities, c_mu, nodal_tke, nodal_epsilon, nodal_f_mu);

        rValues.clear();
        rValues.push_back(nodal_nu_t_sensitivities);
    };

    auto calculate_primal_quantities = [](std::vector<double>& rSensitivities,
                                          const NodeType& rNode,
                                          const ProcessInfo& rCurrentProcessInfo) {
        rSensitivities.push_back(rNode.FastGetSolutionStepValue(TURBULENT_VISCOSITY));
    };

    RansModellingApplicationTestUtilities::RunNodalScalarSensitivityTest(
        r_model_part, y_plus_model_process, calculate_primal_quantities,
        calculate_sensitivities, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        perturb_variable, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonNodalEpsilonSensitivities, RANSEvModelsKEpsilonNodalMatrices)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("RansSensitivities");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelProcess y_plus_model_process(r_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivities = [](std::vector<Vector>& rValues,
                                      const ElementType& rElement,
                                      const ProcessInfo& rCurrentProcessInfo) {
        RansCalculationUtilities rans_calculation_utilities;

        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        const GeometryType& r_geometry = rElement.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Vector nodal_y_plus(number_of_nodes);
        Vector nodal_tke(number_of_nodes);
        Vector nodal_epsilon(number_of_nodes);
        Vector nodal_nu_t(number_of_nodes);
        Vector nodal_f_mu(number_of_nodes);

        RansEvmKEpsilonModel::ReadNodalDataFromElement(
            nodal_y_plus, nodal_tke, nodal_epsilon, nodal_nu_t, nodal_f_mu, rElement);

        Vector nodal_nu_t_sensitivities;
        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
            nodal_nu_t_sensitivities, c_mu, nodal_tke, nodal_epsilon, nodal_f_mu);

        rValues.clear();
        rValues.push_back(nodal_nu_t_sensitivities);
    };

    auto calculate_primal_quantities = [](std::vector<double>& rSensitivities,
                                          const NodeType& rNode,
                                          const ProcessInfo& rCurrentProcessInfo) {
        rSensitivities.push_back(rNode.FastGetSolutionStepValue(TURBULENT_VISCOSITY));
    };

    RansModellingApplicationTestUtilities::RunNodalScalarSensitivityTest(
        r_model_part, y_plus_model_process, calculate_primal_quantities,
        calculate_sensitivities, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        perturb_variable, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonNodalVelocitySensitivities, RANSEvModelsKEpsilonNodalMatrices)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("RansSensitivities");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_model_sensitivities_process(
        r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess y_plus_model_process(r_model_part, empty_parameters);

    y_plus_model_process.Check();
    y_plus_model_process.Execute();

    y_plus_model_sensitivities_process.Check();
    y_plus_model_sensitivities_process.Execute();

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivities = [](std::vector<Matrix>& rValues,
                                      const ElementType& rElement,
                                      const ProcessInfo& rCurrentProcessInfo) {
        RansCalculationUtilities rans_calculation_utilities;

        const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

        const GeometryType& r_geometry = rElement.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        Vector nodal_y_plus(number_of_nodes);
        Vector nodal_tke(number_of_nodes);
        Vector nodal_epsilon(number_of_nodes);
        Vector nodal_nu_t(number_of_nodes);
        Vector nodal_f_mu(number_of_nodes);

        RansEvmKEpsilonModel::ReadNodalDataFromElement(
            nodal_y_plus, nodal_tke, nodal_epsilon, nodal_nu_t, nodal_f_mu, rElement);

        const Matrix& r_y_plus_sensitivities =
            rElement.GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

        Matrix f_mu_velocity_sensitivities;
        EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
            f_mu_velocity_sensitivities, nodal_y_plus, r_y_plus_sensitivities);

        Matrix nodal_nu_t_sensitivities;
        EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
            nodal_nu_t_sensitivities, c_mu, nodal_tke, nodal_epsilon, f_mu_velocity_sensitivities);

        rValues.clear();
        rValues.push_back(nodal_nu_t_sensitivities);
    };

    auto calculate_primal_quantities = [](std::vector<double>& rSensitivities,
                                          const NodeType& rNode,
                                          const ProcessInfo& rCurrentProcessInfo) {
        rSensitivities.push_back(rNode.FastGetSolutionStepValue(TURBULENT_VISCOSITY));
    };

    RansModellingApplicationTestUtilities::RunNodalVectorSensitivityTest(
        r_model_part, y_plus_model_process, calculate_primal_quantities,
        calculate_sensitivities, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        perturb_variable, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementVelocityDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        return r_velocity[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementShapeSensitivity,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode, const int Dim) -> double& {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        return r_coordinates[Dim];
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualVectorSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementTKEFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");

    r_primal_model_part.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 102.3);
    r_adjoint_model_part.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 102.3);

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVMSAdjointElementEpsilonFirstDerivativeLHSMatrix,
                          RANSEvModelsKEpsilonElementResidualMatrices)
{
    Model primal_model;
    ModelPart& r_primal_model_part =
        primal_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_primal_model_part, "VMS2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part =
        adjoint_model.CreateModelPart("RansEvmKElementSensitivityMatrix");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(
        r_adjoint_model_part, "RANSEVMKEpsilonVMSAdjoint2D3N");


    r_primal_model_part.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 102.3);
    r_adjoint_model_part.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 102.3);

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess y_plus_sensitivities_process(
        r_adjoint_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess adjoint_y_plus_process(
        r_adjoint_model_part, empty_parameters);

    RansLogarithmicYPlusModelProcess primal_y_plus_process(r_primal_model_part, empty_parameters);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunElementResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_y_plus_process,
        adjoint_y_plus_process, y_plus_sensitivities_process,
        RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-5, 1e-5);
}
} // namespace Testing
} // namespace Kratos.
