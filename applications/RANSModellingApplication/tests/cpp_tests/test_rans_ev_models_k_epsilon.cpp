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

// Application includes

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
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
}

void InitializeProcessInfo(ModelPart& rModelPart)
{
    // Process info creation
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 0.01);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 0.09);
    rModelPart.GetProcessInfo().SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 1.03);
    rModelPart.GetProcessInfo().SetValue(BOSSAK_ALPHA, -0.03);
    rModelPart.GetProcessInfo().SetValue(NEWMARK_GAMMA, 1.25);

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
    // Set the VELOCITY and PRESSURE nodal values
    array_1d<double, 3> v_1 = ZeroVector(3);
    array_1d<double, 3> v_2 = ZeroVector(3);
    array_1d<double, 3> v_3 = ZeroVector(3);
    v_1[0] = 10.0;
    v_1[1] = 20.0;
    v_2[0] = 15.0;
    v_2[1] = 18.0;
    v_3[0] = 12.0;
    v_3[1] = 11.0;
    (rModelPart.GetNode(1)).GetSolutionStepValue(VELOCITY) = v_1;
    (rModelPart.GetNode(2)).GetSolutionStepValue(VELOCITY) = v_2;
    (rModelPart.GetNode(3)).GetSolutionStepValue(VELOCITY) = v_3;

    (rModelPart.GetNode(1)).GetSolutionStepValue(DISTANCE) = 0.01;
    (rModelPart.GetNode(2)).GetSolutionStepValue(DISTANCE) = 0.02;
    (rModelPart.GetNode(3)).GetSolutionStepValue(DISTANCE) = 0.05;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 10.145;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 12.154;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 18.326;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 5.12;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 1.26;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 6.23;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE, 1) = 4.52;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE, 1) = 8.36;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE, 1) = 2.13;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 9.02;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 6.16;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 4.53;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 6.42;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 3.36;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 4.23;

    // Set the DENSITY and DYNAMIC_VISCOSITY nodal values
    for (ModelPart::NodeIterator it_node = rModelPart.NodesBegin();
         it_node < rModelPart.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 3e-2;
    }
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

        double& nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
        nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke, epsilon, f_mu);

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
    }
}

void CalculateRansYPlusAndUpdateVariables(Process& rProcess, ModelPart& rModelPart)
{
    rProcess.Execute();
    UpdateVariablesInModelPart(rModelPart);
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

void CheckGaussScalarSensitivities(Process& rPrimalProcess,
                                   ModelPart& rModelPart,
                                   Element& rElement,
                                   const Variable<double>& rVariable,
                                   const Vector& rGaussShapeFunctions,
                                   const Matrix& rGaussShapeDerivatives,
                                   const Vector& nu_t_sensitivities,
                                   const Vector& p_k_sensitivities,
                                   const Vector& theta_sensitivities,
                                   const Vector& Re_t_sensitivities,
                                   const Vector& f2_sensitivities)
{
    auto& r_geometry = rElement.GetGeometry();
    int number_of_nodes = static_cast<int>(r_geometry.PointsNumber());

    RansCalculationUtilities rans_calculation_utilities;

    const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    const double nu_t_0 = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_VISCOSITY, rGaussShapeFunctions);
    const double tke_0 = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_KINETIC_ENERGY, rGaussShapeFunctions);
    const double y_plus_0 = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, RANS_Y_PLUS, rGaussShapeFunctions);
    const double epsilon_0 = rans_calculation_utilities.EvaluateInPoint(
        r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, rGaussShapeFunctions);
    const double f_mu_0 = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0);

    BoundedMatrix<double, 2, 2> velocity_gradient_matrix_0;
    rans_calculation_utilities.CalculateGradient<2>(
        velocity_gradient_matrix_0, r_geometry, VELOCITY, rGaussShapeDerivatives);
    const double P_k_0 = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
        velocity_gradient_matrix_0, nu_t_0, tke_0);
    const double theta_0 =
        EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu_0, tke_0, nu_t_0);
    const double Re_t_0 = std::pow(tke_0, 2) / (epsilon_0 * nu_t_0);
    const double f2_0 = EvmKepsilonModelUtilities::CalculateF2(tke_0, nu_t_0, epsilon_0);

    const double delta = 1e-8;
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        double& nodal_value = r_geometry[i_node].FastGetSolutionStepValue(rVariable);
        nodal_value += delta;
        CalculateRansYPlusAndUpdateVariables(rPrimalProcess, rModelPart);

        const double y_plus = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, rGaussShapeFunctions);
        const double nu_t = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_VISCOSITY, rGaussShapeFunctions);
        const double tke = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, rGaussShapeFunctions);
        const double epsilon = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, rGaussShapeFunctions);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix, r_geometry, VELOCITY, rGaussShapeDerivatives);
        const double P_k = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
            velocity_gradient_matrix, nu_t, tke);

        const double theta =
            EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
        const double Re_t = std::pow(tke, 2) / (nu_t * epsilon);
        const double f2 = EvmKepsilonModelUtilities::CalculateF2(tke, nu_t, epsilon);

        const double nu_t_sensitivity = (nu_t - nu_t_0) / delta;
        const double P_k_sensitivity = (P_k - P_k_0) / delta;
        const double theta_sensitivity = (theta - theta_0) / delta;
        const double Re_t_sensitivity = (Re_t - Re_t_0) / delta;
        const double f2_sensitivity = (f2 - f2_0) / delta;

        KRATOS_CHECK_NEAR(nu_t_sensitivities[i_node], nu_t_sensitivity, 1e-6);
        KRATOS_CHECK_NEAR(p_k_sensitivities[i_node], P_k_sensitivity, 1e-5);
        KRATOS_CHECK_NEAR(theta_sensitivities[i_node], theta_sensitivity, 1e-6);
        KRATOS_CHECK_NEAR(Re_t_sensitivities[i_node], Re_t_sensitivity, 1e-5);
        KRATOS_CHECK_NEAR(f2_sensitivities[i_node], f2_sensitivity, 1e-5);

        KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);
        KRATOS_CHECK_NOT_EQUAL(P_k_sensitivity, 0.0);
        KRATOS_CHECK_NOT_EQUAL(theta_sensitivity, 0.0);
        KRATOS_CHECK_NOT_EQUAL(Re_t_sensitivity, 0.0);

        nodal_value -= delta;
    }
}

} // namespace RansEvmKEpsilonModel

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonNodalVelocitySensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonVelocitySensitivitiesNodal");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();
    Matrix& r_adjoint_values = r_element.GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Matrix f_mu_velocity_sensitivities;
    Matrix nu_t_velocity_sensitivities;

    Vector y_plus_0(number_of_nodes);
    Vector tke_0(number_of_nodes);
    Vector epsilon_0(number_of_nodes);
    Vector nu_t_0(number_of_nodes);
    Vector f_mu_0(number_of_nodes);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(y_plus_0, tke_0, epsilon_0,
                                                   nu_t_0, f_mu_0, r_element);

    EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
        f_mu_velocity_sensitivities, y_plus_0, r_adjoint_values);

    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
        nu_t_velocity_sensitivities, c_mu, tke_0, epsilon_0, f_mu_velocity_sensitivities);

    // Calculate finite difference values
    const double delta = 1e-8;
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        array_1d<double, 3>& r_velocity =
            r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
        for (int i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            r_velocity[i_dim] += delta;
            RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(
                primal_process, r_model_part);

            const double nu_t =
                r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            const double nu_t_sensitivity = (nu_t - nu_t_0[i_node]) / delta;
            KRATOS_CHECK_NEAR(nu_t_velocity_sensitivities(i_node, i_dim),
                              nu_t_sensitivity, 1e-6);
            KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);

            r_velocity[i_dim] -= delta;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonNodalTKESensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonTKESensitivitiesNodal");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();

    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector f_mu_0(number_of_nodes);
    Vector y_plus_0(number_of_nodes);
    Vector tke_0(number_of_nodes);
    Vector epsilon_0(number_of_nodes);
    Vector nu_t_0(number_of_nodes);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(y_plus_0, tke_0, epsilon_0,
                                                   nu_t_0, f_mu_0, r_element);

    Vector nu_t_velocity_sensitivities;
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
        nu_t_velocity_sensitivities, c_mu, tke_0, epsilon_0, f_mu_0);

    // Calculate finite difference values
    const double delta = 1e-8;
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        double& tke_nodal =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        tke_nodal += delta;
        RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

        const double nu_t = r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        const double nu_t_sensitivity = (nu_t - nu_t_0[i_node]) / delta;

        KRATOS_CHECK_NEAR(nu_t_velocity_sensitivities[i_node], nu_t_sensitivity, 1e-6);
        KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);

        tke_nodal -= delta;
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonNodalEpsilonSensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonEpsilonSensitivitiesNodal");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();

    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector f_mu_0(number_of_nodes);
    Vector y_plus_0(number_of_nodes);
    Vector tke_0(number_of_nodes);
    Vector epsilon_0(number_of_nodes);
    Vector nu_t_0(number_of_nodes);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(y_plus_0, tke_0, epsilon_0,
                                                   nu_t_0, f_mu_0, r_element);

    Vector nu_t_velocity_sensitivities;
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
        nu_t_velocity_sensitivities, c_mu, tke_0, epsilon_0, f_mu_0);

    // Calculate finite difference values
    const double delta = 1e-8;
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        double& epsilon_nodal =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        epsilon_nodal += delta;
        RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

        const double nu_t = r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        const double nu_t_sensitivity = (nu_t - nu_t_0[i_node]) / delta;

        KRATOS_CHECK_NEAR(nu_t_velocity_sensitivities[i_node], nu_t_sensitivity, 1e-6);
        KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);

        epsilon_nodal -= delta;
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussVelocitySensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonVelocitySensitivitiesGauss");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();
    Matrix& r_adjoint_values = r_element.GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    rans_calculation_utilities.CalculateGeometryData(
        r_geometry, r_element.GetIntegrationMethod(), gauss_weights,
        shape_functions, shape_derivatives);

    const int num_gauss_points = gauss_weights.size();

    Vector nodal_y_plus_0(number_of_nodes);
    Vector nodal_tke_0(number_of_nodes);
    Vector nodal_epsilon_0(number_of_nodes);
    Vector nodal_nu_t_0(number_of_nodes);
    Vector nodal_f_mu_0(number_of_nodes);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(
        nodal_y_plus_0, nodal_tke_0, nodal_epsilon_0, nodal_nu_t_0, nodal_f_mu_0, r_element);

    Matrix nodal_f_mu_sensitivities(number_of_nodes, domain_size);
    EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVectorSensitivities(
        nodal_f_mu_sensitivities, nodal_y_plus_0, r_adjoint_values);

    Matrix nodal_nu_t_sensitivities(number_of_nodes, domain_size);
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVectorSensitivities(
        nodal_nu_t_sensitivities, c_mu, nodal_tke_0, nodal_epsilon_0, nodal_f_mu_sensitivities);

    // Calculate finite difference values
    const double delta = 1e-8;
    for (int g = 1; g < num_gauss_points; ++g)
    {
        RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const Matrix& r_shape_derivatives = shape_derivatives[g];

        const double y_plus_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, gauss_shape_functions);
        const double tke_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double nu_t_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix_0;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix_0, r_geometry, VELOCITY, r_shape_derivatives);

        const double f_mu_0 = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0);
        const double P_k_0 = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
            velocity_gradient_matrix_0, nu_t_0, tke_0);
        const double theta_0 =
            EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu_0, tke_0, nu_t_0);

        const double Re_t_0 = std::pow(tke_0, 2) / (nu_t_0 * epsilon_0);
        const double f2_0 = EvmKepsilonModelUtilities::CalculateF2(tke_0, nu_t_0, epsilon_0);

        Matrix gauss_f_mu_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussFmuVectorSensitivities(
            gauss_f_mu_velocity_sensitivities, y_plus_0, r_adjoint_values, gauss_shape_functions);

        Matrix gauss_nu_t_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_nu_t_velocity_sensitivities, nodal_nu_t_sensitivities, gauss_shape_functions);

        Matrix gauss_y_plus_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_y_plus_velocity_sensitivities, r_adjoint_values, gauss_shape_functions);

        Matrix gauss_production_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateProductionVelocitySensitivities<2>(
            gauss_production_velocity_sensitivities, nu_t_0, gauss_nu_t_velocity_sensitivities,
            velocity_gradient_matrix_0, r_shape_derivatives);

        Matrix gauss_theta_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateThetaVelocitySensitivity(
            gauss_theta_sensitivities, c_mu, f_mu_0, tke_0, nu_t_0,
            gauss_f_mu_velocity_sensitivities, gauss_nu_t_velocity_sensitivities);

        Matrix gauss_turbulent_reynolds_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateTurbulentReynoldsNumberVelocitySensitivity(
            gauss_turbulent_reynolds_sensitivities, tke_0, epsilon_0, nu_t_0,
            gauss_nu_t_velocity_sensitivities);

        Matrix gauss_f2_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateF2VelocitySensitivity(
            gauss_f2_sensitivities, tke_0, epsilon_0, nu_t_0, gauss_nu_t_velocity_sensitivities);

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            array_1d<double, 3>& r_velocity =
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
            for (int i_dim = 0; i_dim < domain_size; ++i_dim)
            {
                r_velocity[i_dim] += delta;
                RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(
                    primal_process, r_model_part);

                const double y_plus = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, RANS_Y_PLUS, gauss_shape_functions);
                const double nu_t = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);
                const double tke = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
                const double epsilon = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);

                const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);

                BoundedMatrix<double, 2, 2> velocity_gradient_matrix;
                rans_calculation_utilities.CalculateGradient<2>(
                    velocity_gradient_matrix, r_geometry, VELOCITY, r_shape_derivatives);

                const double P_k = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
                    velocity_gradient_matrix, nu_t, tke);

                const double theta =
                    EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
                const double Re_t = std::pow(tke, 2) / (nu_t * epsilon);
                const double f2 =
                    EvmKepsilonModelUtilities::CalculateF2(tke, nu_t, epsilon);

                const double f_mu_sensitivity = (f_mu - f_mu_0) / delta;
                const double nu_t_sensitivity = (nu_t - nu_t_0) / delta;
                const double y_plus_sensitivity = (y_plus - y_plus_0) / delta;
                const double P_k_sensitivity = (P_k - P_k_0) / delta;
                const double theta_sensitivity = (theta - theta_0) / delta;
                const double Re_t_sensitivity = (Re_t - Re_t_0) / delta;
                const double f2_sensitivity = (f2 - f2_0) / delta;

                KRATOS_CHECK_NEAR(gauss_f_mu_velocity_sensitivities(i_node, i_dim),
                                  f_mu_sensitivity, 1e-6);
                KRATOS_CHECK_NEAR(gauss_y_plus_velocity_sensitivities(i_node, i_dim),
                                  y_plus_sensitivity, 1e-6);
                KRATOS_CHECK_NEAR(gauss_nu_t_velocity_sensitivities(i_node, i_dim),
                                  nu_t_sensitivity, 1e-6);
                KRATOS_CHECK_NEAR(gauss_production_velocity_sensitivities(i_node, i_dim),
                                  P_k_sensitivity, 1e-5);
                KRATOS_CHECK_NEAR(gauss_theta_sensitivities(i_node, i_dim),
                                  theta_sensitivity, 1e-6);
                KRATOS_CHECK_NEAR(gauss_turbulent_reynolds_sensitivities(i_node, i_dim),
                                  Re_t_sensitivity, 1e-5);
                KRATOS_CHECK_NEAR(gauss_f2_sensitivities(i_node, i_dim),
                                  f2_sensitivity, 1e-6);

                KRATOS_CHECK_NOT_EQUAL(f_mu_sensitivity, 0.0);
                KRATOS_CHECK_NOT_EQUAL(y_plus_sensitivity, 0.0);
                KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);
                KRATOS_CHECK_NOT_EQUAL(P_k_sensitivity, 0.0);
                KRATOS_CHECK_NOT_EQUAL(theta_sensitivity, 0.0);
                KRATOS_CHECK_NOT_EQUAL(Re_t_sensitivity, 0.0);

                r_velocity[i_dim] -= delta;
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussTKESensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonTKESensitivitiesGauss");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    rans_calculation_utilities.CalculateGeometryData(
        r_geometry, r_element.GetIntegrationMethod(), gauss_weights,
        shape_functions, shape_derivatives);

    const int num_gauss_points = gauss_weights.size();

    Vector nodal_y_plus_0(number_of_nodes);
    Vector nodal_tke_0(number_of_nodes);
    Vector nodal_epsilon_0(number_of_nodes);
    Vector nodal_nu_t_0(number_of_nodes);
    Vector nodal_f_mu_0(number_of_nodes);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(
        nodal_y_plus_0, nodal_tke_0, nodal_epsilon_0, nodal_nu_t_0, nodal_f_mu_0, r_element);

    Vector nodal_nu_t_sensitivities(number_of_nodes, domain_size);
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
        nodal_nu_t_sensitivities, c_mu, nodal_tke_0, nodal_epsilon_0, nodal_f_mu_0);

    // Calculate finite difference values
    for (int g = 0; g < num_gauss_points; ++g)
    {
        RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const Matrix& r_shape_derivatives = shape_derivatives[g];

        const double y_plus_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, gauss_shape_functions);
        const double tke_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double nu_t_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix_0;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix_0, r_geometry, VELOCITY, r_shape_derivatives);

        const double f_mu_0 = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0);
        const double Re_t_0 = std::pow(tke_0, 2) / (nu_t_0 * epsilon_0);

        Vector gauss_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_nu_t_sensitivities, nodal_nu_t_sensitivities, gauss_shape_functions);

        Vector gauss_production_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<2>(
            gauss_production_sensitivities, gauss_nu_t_sensitivities, velocity_gradient_matrix_0);

        Vector gauss_theta_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateThetaTKESensitivity(
            gauss_theta_sensitivities, c_mu, f_mu_0, tke_0, nu_t_0,
            gauss_nu_t_sensitivities, gauss_shape_functions);

        Vector gauss_re_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateTurbulentReynoldsNumberTKESensitivity(
            gauss_re_t_sensitivities, tke_0, epsilon_0, nu_t_0,
            gauss_nu_t_sensitivities, gauss_shape_functions);

        Vector gauss_f2_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateF2ScalarSensitivity(
            gauss_f2_sensitivities, epsilon_0, Re_t_0, gauss_re_t_sensitivities);

        RansEvmKEpsilonModel::CheckGaussScalarSensitivities(
            primal_process, r_model_part, r_element, TURBULENT_KINETIC_ENERGY,
            gauss_shape_functions, r_shape_derivatives, gauss_nu_t_sensitivities,
            gauss_production_sensitivities, gauss_theta_sensitivities,
            gauss_re_t_sensitivities, gauss_f2_sensitivities);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussEpsilonSensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonEpsilonSensitivitiesGauss");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    rans_calculation_utilities.CalculateGeometryData(
        r_geometry, r_element.GetIntegrationMethod(), gauss_weights,
        shape_functions, shape_derivatives);

    const int num_gauss_points = gauss_weights.size();

    Vector nodal_y_plus_0(number_of_nodes);
    Vector nodal_tke_0(number_of_nodes);
    Vector nodal_epsilon_0(number_of_nodes);
    Vector nodal_nu_t_0(number_of_nodes);
    Vector nodal_f_mu_0(number_of_nodes);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(
        nodal_y_plus_0, nodal_tke_0, nodal_epsilon_0, nodal_nu_t_0, nodal_f_mu_0, r_element);

    Vector nodal_nu_t_sensitivities(number_of_nodes, domain_size);
    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
        nodal_nu_t_sensitivities, c_mu, nodal_tke_0, nodal_epsilon_0, nodal_f_mu_0);

    // Calculate finite difference values
    for (int g = 0; g < num_gauss_points; ++g)
    {
        RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);
        const Vector& gauss_shape_functions = row(shape_functions, g);
        const Matrix& r_shape_derivatives = shape_derivatives[g];

        const double y_plus_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, gauss_shape_functions);
        const double tke_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);
        const double nu_t_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);
        const double f_mu_0 = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0);
        const double Re_t_0 = std::pow(tke_0, 2) / (nu_t_0 * epsilon_0);

        Vector gauss_nu_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateGaussSensitivities(
            gauss_nu_t_sensitivities, nodal_nu_t_sensitivities, gauss_shape_functions);

        BoundedMatrix<double, 2, 2> velocity_gradient_matrix_0;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_matrix_0, r_geometry, VELOCITY, r_shape_derivatives);

        Vector gauss_production_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateProductionScalarSensitivities<2>(
            gauss_production_sensitivities, gauss_nu_t_sensitivities, velocity_gradient_matrix_0);

        Vector gauss_theta_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateThetaEpsilonSensitivity(
            gauss_theta_sensitivities, c_mu, f_mu_0, tke_0, nu_t_0, gauss_nu_t_sensitivities);

        Vector gauss_re_t_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateTurbulentReynoldsNumberEpsilonSensitivity(
            gauss_re_t_sensitivities, tke_0, epsilon_0, nu_t_0,
            gauss_nu_t_sensitivities, gauss_shape_functions);

        Vector gauss_f2_sensitivities(number_of_nodes);
        EvmKepsilonModelAdjointUtilities::CalculateF2ScalarSensitivity(
            gauss_f2_sensitivities, epsilon_0, Re_t_0, gauss_re_t_sensitivities);

        RansEvmKEpsilonModel::CheckGaussScalarSensitivities(
            primal_process, r_model_part, r_element, TURBULENT_ENERGY_DISSIPATION_RATE,
            gauss_shape_functions, r_shape_derivatives, gauss_nu_t_sensitivities,
            gauss_production_sensitivities, gauss_theta_sensitivities,
            gauss_re_t_sensitivities, gauss_f2_sensitivities);

        // checking pde terms
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonGaussShapeSensitivities, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonShapeSensitivitiesGauss");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonTestModelPart(r_model_part,
                                                               "RANSEVMK2D3N");

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector gauss_weights_0;
    Matrix shape_functions_0;
    ShapeFunctionDerivativesArrayType shape_derivatives_0;
    rans_calculation_utilities.CalculateGeometryData(
        r_geometry, r_element.GetIntegrationMethod(), gauss_weights_0,
        shape_functions_0, shape_derivatives_0);

    const int num_gauss_points = gauss_weights_0.size();

    Vector nodal_y_plus_0(number_of_nodes);
    Vector nodal_tke_0(number_of_nodes);
    Vector nodal_epsilon_0(number_of_nodes);
    Vector nodal_nu_t_0(number_of_nodes);
    Vector nodal_f_mu_0(number_of_nodes);
    Matrix nodal_velocity(number_of_nodes, domain_size);

    RansEvmKEpsilonModel::ReadNodalDataFromElement(
        nodal_y_plus_0, nodal_tke_0, nodal_epsilon_0, nodal_nu_t_0, nodal_f_mu_0, r_element);

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        for (int i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            nodal_velocity(i_node, i_dim) =
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY)[i_dim];
        }
    }

    Geometry<Point>::JacobiansType J;
    r_geometry.Jacobian(J, r_element.GetIntegrationMethod());

    Geometry<Point>::ShapeFunctionsGradientsType DN_De;
    DN_De = r_geometry.ShapeFunctionsLocalGradients(r_element.GetIntegrationMethod());

    // Calculate finite difference values
    const double delta = 1e-8;
    for (int g = 1; g < num_gauss_points; ++g)
    {
        RansEvmKEpsilonModel::CalculateRansYPlusAndUpdateVariables(primal_process, r_model_part);
        const Vector& gauss_shape_functions_0 = row(shape_functions_0, g);
        const Matrix& r_shape_derivatives_0 = shape_derivatives_0[g];

        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        ShapeParameter deriv;

        const double nu_t_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions_0);
        const double tke_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions_0);

        BoundedMatrix<double, 2, 2> velocity_gradient_0;
        rans_calculation_utilities.CalculateGradient<2>(
            velocity_gradient_0, r_geometry, VELOCITY, r_shape_derivatives_0);

        // calculating initial fd values
        const double production_0 = EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
            velocity_gradient_0, nu_t_0, tke_0);

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            for (int i_dim = 0; i_dim < domain_size; ++i_dim)
            {
                deriv.NodeIndex = i_node;
                deriv.Direction = i_dim;

                double detJ_deriv;
                GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;
                geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);

                double nu_t_deriv = 0.0;

                double production_adjoint_derivative;
                EvmKepsilonModelAdjointUtilities::CalculateProductionShapeSensitivities<2>(
                    production_adjoint_derivative, nu_t_0, nu_t_deriv, nodal_velocity,
                    velocity_gradient_0, r_shape_derivatives_0, DN_DX_deriv);

                // doing fd sensitivities
                if (i_dim == 0)
                    r_geometry[i_node].X() += delta;
                else if (i_dim == 1)
                    r_geometry[i_node].Y() += delta;

                Vector gauss_weights;
                Matrix shape_functions;
                ShapeFunctionDerivativesArrayType shape_derivatives;
                rans_calculation_utilities.CalculateGeometryData(
                    r_geometry, r_element.GetIntegrationMethod(), gauss_weights,
                    shape_functions, shape_derivatives);

                const Vector& gauss_shape_functions = row(shape_functions, g);
                const Matrix& r_shape_derivatives = shape_derivatives[g];

                BoundedMatrix<double, 2, 2> velocity_gradient;
                rans_calculation_utilities.CalculateGradient<2>(
                    velocity_gradient, r_geometry, VELOCITY, r_shape_derivatives);

                const double nu_t = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, TURBULENT_VISCOSITY, gauss_shape_functions);
                const double tke = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);

                const double production =
                    EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
                        velocity_gradient, nu_t, tke);

                const double production_sensitivity = (production - production_0) / delta;

                KRATOS_CHECK_NEAR(production_sensitivity,
                                  production_adjoint_derivative, 1e-5);
                KRATOS_CHECK_NOT_EQUAL(production_sensitivity, 0.0);

                if (i_dim == 0)
                    r_geometry[i_node].X() -= delta;
                else if (i_dim == 1)
                    r_geometry[i_node].Y() -= delta;
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementSensitivityMatrix, RANSEvModels)
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

    auto perturb_variable = [](NodeType& rNode, const int iDim, const double Epsilon) {
        array_1d<double, 3>& r_coordinates = rNode.Coordinates();
        r_coordinates[iDim] += Epsilon;
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        TURBULENT_KINETIC_ENERGY, RANS_AUXILIARY_VARIABLE_1, r_primal_model_part,
        r_adjoint_model_part, primal_y_plus_process, adjoint_y_plus_process,
        y_plus_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementVelocityFirstDerivativeLHSMatrix, RANSEvModels)
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

    auto perturb_variable = [](NodeType& rNode, const int iDim, const double Epsilon) {
        array_1d<double, 3>& nodal_value = rNode.FastGetSolutionStepValue(VELOCITY);
        nodal_value[iDim] += Epsilon;
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           const ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_VELOCITY_PARTIAL_DERIVATIVE, rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualVectorSensitivityTest(
        TURBULENT_KINETIC_ENERGY, RANS_AUXILIARY_VARIABLE_1, r_primal_model_part,
        r_adjoint_model_part, primal_y_plus_process, adjoint_y_plus_process,
        y_plus_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementTKEFirstDerivativeLHSMatrix, RANSEvModels)
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

    auto perturb_variable = [](NodeType& rNode, const double Epsilon) {
        double& value = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        value += Epsilon;
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        TURBULENT_KINETIC_ENERGY, RANS_AUXILIARY_VARIABLE_1, r_primal_model_part,
        r_adjoint_model_part, primal_y_plus_process, adjoint_y_plus_process,
        y_plus_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementEpsilonFirstDerivativeLHSMatrix, RANSEvModels)
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

    auto perturb_variable = [](NodeType& rNode, const double Epsilon) {
        double& value = rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        value += Epsilon;
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.Calculate(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE,
                           rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        TURBULENT_KINETIC_ENERGY, RANS_AUXILIARY_VARIABLE_1, r_primal_model_part,
        r_adjoint_model_part, primal_y_plus_process, adjoint_y_plus_process,
        y_plus_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKElementTKESecondDerivativeLHSMatrix, RANSEvModels)
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

    auto perturb_variable = [](NodeType& rNode, const double Epsilon) {
        double& value = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE);
        value += Epsilon;
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rOutput, rProcessInfo);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        TURBULENT_KINETIC_ENERGY, RANS_AUXILIARY_VARIABLE_1, r_primal_model_part,
        r_adjoint_model_part, primal_y_plus_process, adjoint_y_plus_process,
        y_plus_sensitivities_process, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos.
