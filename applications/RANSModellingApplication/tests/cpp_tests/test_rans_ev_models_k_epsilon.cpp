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
#include "includes/model_part.h"
#include "rans_modelling_application_variables.h"
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
void GenerateRansEvmKEpsilonTestModelPart(ModelPart& rModelPart)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);

    // Process info creation
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(TURBULENCE_RANS_C_MU, 0.09);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 3.0e-02);

    // Element creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
    rModelPart.CreateNewElement("RANSEVMK2D3N", 1, elem_nodes, p_elem_prop);

    // Set the VELOCITY and PRESSURE nodal values
    array_1d<double, 3> v_1 = ZeroVector(3);
    array_1d<double, 3> v_2 = ZeroVector(3);
    array_1d<double, 3> v_3 = ZeroVector(3);
    v_1[0] = 10.0;
    v_1[1] = 20.0;
    v_2[0] = 1000.0;
    v_2[1] = 500.0;
    v_3[0] = 30.0;
    v_3[1] = 20.0;
    (rModelPart.GetNode(1)).GetSolutionStepValue(VELOCITY) = v_1;
    (rModelPart.GetNode(2)).GetSolutionStepValue(VELOCITY) = v_2;
    (rModelPart.GetNode(3)).GetSolutionStepValue(VELOCITY) = v_3;

    (rModelPart.GetNode(1)).GetSolutionStepValue(DISTANCE) = 0.01;
    (rModelPart.GetNode(2)).GetSolutionStepValue(DISTANCE) = 0.02;
    (rModelPart.GetNode(3)).GetSolutionStepValue(DISTANCE) = 0.05;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 0.145;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 200.154;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 20.326;

    (rModelPart.GetNode(1)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 5.124;
    (rModelPart.GetNode(2)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 150.26;
    (rModelPart.GetNode(3)).GetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 54.23;

    // Set the DENSITY and DYNAMIC_VISCOSITY nodal values
    for (ModelPart::NodeIterator it_node = rModelPart.NodesBegin();
         it_node < rModelPart.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY) =
            p_elem_prop->GetValue(KINEMATIC_VISCOSITY);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVelocitySensitivitiesNodal, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonVelocitySensitivitiesNodal");
    GenerateRansEvmKEpsilonTestModelPart(r_model_part);

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate finite difference values

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    primal_process.Execute();

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();
    Matrix& r_adjoint_values = r_element.GetValue(RANS_Y_PLUS_SENSITIVITIES);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Matrix f_mu_velocity_sensitivities;
    Matrix nu_t_velocity_sensitivities;

    const double delta = 1e-8;
    primal_process.Execute();

    Vector f_mu_0(number_of_nodes);
    Vector y_plus_0(number_of_nodes);
    Vector tke_0(number_of_nodes);
    Vector epsilon_0(number_of_nodes);
    Vector nu_t_0(number_of_nodes);
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        tke_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        epsilon_0[i_node] =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        y_plus_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
        f_mu_0[i_node] = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0[i_node]);
        nu_t_0[i_node] = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke_0[i_node], epsilon_0[i_node], f_mu_0[i_node]);
    }

    EvmKepsilonModelAdjointUtilities::CalculateNodalFmuVelocitySensitivities(
        f_mu_velocity_sensitivities, y_plus_0, r_adjoint_values);

    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityVelocitySensitivities(
        nu_t_velocity_sensitivities, c_mu, tke_0, epsilon_0, f_mu_velocity_sensitivities);

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        for (int i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            array_1d<double, 3>& r_velocity =
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
            r_velocity[i_dim] += delta;
            primal_process.Execute();

            const double y_plus = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
            const double tke =
                r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double epsilon = r_geometry[i_node].FastGetSolutionStepValue(
                TURBULENT_ENERGY_DISSIPATION_RATE);

            const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
            const double f_mu_sensitivity = (f_mu - f_mu_0[i_node]) / delta;

            const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
                c_mu, tke, epsilon, f_mu);
            const double nu_t_sensitivity = (nu_t - nu_t_0[i_node]) / delta;

            KRATOS_CHECK_NEAR(f_mu_velocity_sensitivities(i_node, i_dim),
                              f_mu_sensitivity, 1e-6);

            KRATOS_CHECK_NEAR(nu_t_velocity_sensitivities(i_node, i_dim),
                              nu_t_sensitivity, 1e-6);

            KRATOS_CHECK_NOT_EQUAL(f_mu_sensitivity, 0.0);
            KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);

            r_velocity[i_dim] -= delta;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonTKESensitivitiesNodal, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonTKESensitivitiesNodal");
    GenerateRansEvmKEpsilonTestModelPart(r_model_part);

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate finite difference values

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    primal_process.Execute();

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();

    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector nu_t_velocity_sensitivities;

    const double delta = 1e-8;
    primal_process.Execute();

    Vector f_mu_0(number_of_nodes);
    Vector y_plus_0(number_of_nodes);
    Vector tke_0(number_of_nodes);
    Vector epsilon_0(number_of_nodes);
    Vector nu_t_0(number_of_nodes);
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        tke_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        epsilon_0[i_node] =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        y_plus_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
        f_mu_0[i_node] = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0[i_node]);
        nu_t_0[i_node] = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke_0[i_node], epsilon_0[i_node], f_mu_0[i_node]);
    }

    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityTKESensitivities(
        nu_t_velocity_sensitivities, c_mu, tke_0, epsilon_0, f_mu_0);

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        double& tke_nodal =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        tke_nodal += delta;
        primal_process.Execute();

        const double y_plus = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
        const double tke =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);

        const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke, epsilon, f_mu);
        const double nu_t_sensitivity = (nu_t - nu_t_0[i_node]) / delta;

        KRATOS_CHECK_NEAR(nu_t_velocity_sensitivities[i_node], nu_t_sensitivity, 1e-6);

        KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);

        tke_nodal -= delta;
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonSensitivitiesNodal, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonEpsilonSensitivitiesNodal");
    GenerateRansEvmKEpsilonTestModelPart(r_model_part);

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate finite difference values

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    primal_process.Execute();

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();

    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector nu_t_velocity_sensitivities;

    const double delta = 1e-8;
    primal_process.Execute();

    Vector f_mu_0(number_of_nodes);
    Vector y_plus_0(number_of_nodes);
    Vector tke_0(number_of_nodes);
    Vector epsilon_0(number_of_nodes);
    Vector nu_t_0(number_of_nodes);
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        tke_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        epsilon_0[i_node] =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        y_plus_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
        f_mu_0[i_node] = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0[i_node]);
        nu_t_0[i_node] = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke_0[i_node], epsilon_0[i_node], f_mu_0[i_node]);
    }

    EvmKepsilonModelAdjointUtilities::CalculateNodalTurbulentViscosityEpsilonSensitivities(
        nu_t_velocity_sensitivities, c_mu, tke_0, epsilon_0, f_mu_0);

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        double& epsilon_nodal =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        epsilon_nodal += delta;
        primal_process.Execute();

        const double y_plus = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
        const double tke =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_geometry[i_node].FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);

        const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            c_mu, tke, epsilon, f_mu);
        const double nu_t_sensitivity = (nu_t - nu_t_0[i_node]) / delta;

        KRATOS_CHECK_NEAR(nu_t_velocity_sensitivities[i_node], nu_t_sensitivity, 1e-6);

        KRATOS_CHECK_NOT_EQUAL(nu_t_sensitivity, 0.0);

        epsilon_nodal -= delta;
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVelocitySensitivitiesGauss, RANSEvModels)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RansEvmKEpsilonVelocitySensitivitiesGauss");
    GenerateRansEvmKEpsilonTestModelPart(r_model_part);

    const double c_mu = r_model_part.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    Parameters empty_parameters = Parameters(R"({})");

    RansLogarithmicYPlusModelSensitivitiesProcess adjoint_process(r_model_part, empty_parameters);
    RansLogarithmicYPlusModelProcess primal_process(r_model_part, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate finite difference values

    // Calculate initial y_plus values
    RansCalculationUtilities rans_calculation_utilities;
    primal_process.Check();
    primal_process.Execute();

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();
    Matrix& r_adjoint_values = r_element.GetValue(RANS_Y_PLUS_SENSITIVITIES);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    rans_calculation_utilities.CalculateGeometryData(
        r_geometry, r_element.GetIntegrationMethod(), gauss_weights,
        shape_functions, shape_derivatives);
    const int num_gauss_points = gauss_weights.size();

    std::vector<Matrix> f_mu_velocity_sensitivities(num_gauss_points);


    const double delta = 1e-8;
    for (int g = 0; g < num_gauss_points; ++g)
    {
        primal_process.Execute();
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const double y_plus_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, RANS_Y_PLUS, gauss_shape_functions);
        const double tke_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_KINETIC_ENERGY, gauss_shape_functions);
        const double epsilon_0 = rans_calculation_utilities.EvaluateInPoint(
            r_geometry, TURBULENT_ENERGY_DISSIPATION_RATE, gauss_shape_functions);

        Matrix gauss_f_mu_velocity_sensitivities(number_of_nodes, domain_size);
        EvmKepsilonModelAdjointUtilities::CalculateGaussFmuVelocitySensitivities(
            gauss_f_mu_velocity_sensitivities, y_plus_0, r_adjoint_values, gauss_shape_functions);

        const double f_mu_0 = EvmKepsilonModelUtilities::CalculateFmu(y_plus_0);

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            for (int i_dim = 0; i_dim < domain_size; ++i_dim)
            {
                array_1d<double, 3>& r_velocity =
                    r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
                r_velocity[i_dim] += delta;
                primal_process.Execute();

                const double y_plus = rans_calculation_utilities.EvaluateInPoint(
                    r_geometry, RANS_Y_PLUS, gauss_shape_functions);

                const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
                const double f_mu_sensitivity = (f_mu - f_mu_0) / delta;


                KRATOS_CHECK_NEAR(gauss_f_mu_velocity_sensitivities(i_node, i_dim),
                                  f_mu_sensitivity, 1e-6);

                KRATOS_CHECK_NOT_EQUAL(f_mu_sensitivity, 0.0);

                r_velocity[i_dim] -= delta;
            }
        }
    }
}
} // namespace Testing
} // namespace Kratos.
