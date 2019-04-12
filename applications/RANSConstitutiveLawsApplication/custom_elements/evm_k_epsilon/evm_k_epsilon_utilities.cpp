#include "evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/cfd_variables.h"
#include "rans_constitutive_laws_application_variables.h"
#include <cmath>
#include <iostream>
#include <limits>

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

double EvmKepsilonModelUtilities::CalculateTurbulentViscosity(const double C_mu,
                                                              const double turbulent_kinetic_energy,
                                                              const double turbulent_energy_dissipation_rate,
                                                              const double f_mu)
{
    return C_mu * f_mu * std::pow(turbulent_kinetic_energy, 2) / turbulent_energy_dissipation_rate;
}

double EvmKepsilonModelUtilities::CalculateFmu(const double y_plus)
{
    return 1.0 - std::exp(-0.0115 * y_plus);
}

double EvmKepsilonModelUtilities::CalculateF2(const double turbulent_kinetic_energy,
                                              const double kinematic_viscosity,
                                              const double turbulent_energy_dissipation_rate)
{
    if (turbulent_energy_dissipation_rate == 0.0)
        return 1.0;
    else
    {
        const double Re_t = std::pow(turbulent_kinetic_energy, 2) /
                            (kinematic_viscosity * turbulent_energy_dissipation_rate);
        const double f2 = 1.0 - 0.22 * std::exp(-1.0 * std::pow(Re_t / 6.0, 2));

        return f2;
    }
}

template <unsigned int TDim>
double EvmKepsilonModelUtilities::CalculateSourceTerm(
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const double turbulent_kinematic_viscosity,
    const double turbulent_kinetic_energy)
{
    const double velocity_divergence =
        RansCalculationUtilities().CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;

    noalias(reynolds_stress_tensor) =
        turbulent_kinematic_viscosity *
        (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity);

    double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            source += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    return source;
}

double EvmKepsilonModelUtilities::CalculateGamma(const double C_mu,
                                                 const double f_mu,
                                                 const double turbulent_kinetic_energy,
                                                 const double turbulent_kinematic_viscosity)
{
    return std::max<double>(
        0.0, C_mu * f_mu * turbulent_kinetic_energy / turbulent_kinematic_viscosity);
}

void EvmKepsilonModelUtilities::CalculatePositiveValuesList(Vector& rOutput, const Vector& rInput)
{
    const int n = rInput.size();
    rOutput.resize(n);

    for (int i = 0; i < n; ++i)
    {
        if (rInput[i] >= std::numeric_limits<double>::epsilon())
            rOutput[i] = 1.0;
        else
            rOutput[i] = 0.0;
    }
}

void EvmKepsilonModelUtilities::CalculateTurbulentValues(double& turbulent_kinetic_energy,
                                                         double& turbulent_energy_dissipation_rate,
                                                         const double y_plus,
                                                         const double kinematic_viscosity,
                                                         const double wall_distance,
                                                         const double c_mu,
                                                         const double von_karman)
{
    const double u_tau = y_plus * kinematic_viscosity / wall_distance;
    turbulent_kinetic_energy = std::pow(u_tau, 2) / std::sqrt(c_mu);
    turbulent_energy_dissipation_rate = std::pow(u_tau, 3) / (von_karman * wall_distance);
}

void EvmKepsilonModelUtilities::CalculateTurbulentValues(double& turbulent_kinetic_energy,
                                                         double& turbulent_energy_dissipation_rate,
                                                         const double velocity_mag,
                                                         const double turbulence_intensity,
                                                         const double mixing_length,
                                                         const double c_mu)
{
    turbulent_kinetic_energy = 1.5 * std::pow(velocity_mag * turbulence_intensity, 2);
    turbulent_energy_dissipation_rate =
        c_mu * std::pow(turbulent_kinetic_energy, 1.5) / mixing_length;
}

// void EvmKepsilonModelUtilities::CalculateTurbulentViscosityForModelPart(ModelPart& rModelPart)
// {
//     KRATOS_TRY

//     int number_of_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
//     const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];

// #pragma omp parallel for
//     for (int i = 0; i < number_of_nodes; ++i)
//     {
//         auto& r_node = *(rModelPart.GetCommunicator().LocalMesh().NodesBegin() + i);

//         const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
//         const double epsilon =
//             r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
//         const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);

//         const double f_mu = CalculateFmu(y_plus);

//         const double nu_t = CalculateTurbulentViscosity(c_mu, tke, epsilon,
//         f_mu); r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
//     }

//     KRATOS_CATCH("");
// }

void EvmKepsilonModelUtilities::UpdateBoundaryConditions(ModelPart& rModelPart)
{
    KRATOS_TRY

    int number_of_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];
    const double von_karman = r_process_info[WALL_VON_KARMAN];

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; ++i)
    {
        auto& r_node = *(rModelPart.GetCommunicator().LocalMesh().NodesBegin() + i);

        if (r_node.Is(INLET))
        {
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);
            double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            double& epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

            CalculateTurbulentValues(tke, epsilon, y_plus, nu, wall_distance, c_mu, von_karman);
        }
        else if (r_node.Is(STRUCTURE))
        {
            double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            double& epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

            tke = 1e-10;
            epsilon = 1e-10;
        }
    }

    KRATOS_CATCH("");
}

void EvmKepsilonModelUtilities::AssignInitialValues(ModelPart& rModelPart)
{
    KRATOS_TRY

    int number_of_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];
    const double von_karman = r_process_info[WALL_VON_KARMAN];

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; ++i)
    {
        auto& r_node = *(rModelPart.GetCommunicator().LocalMesh().NodesBegin() + i);

        if (!r_node.IsFixed(TURBULENT_KINETIC_ENERGY) &&
            !r_node.IsFixed(TURBULENT_ENERGY_DISSIPATION_RATE))
        {
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);
            double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            double& epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

            CalculateTurbulentValues(tke, epsilon, y_plus, nu, wall_distance, c_mu, von_karman);
        }
    }

    KRATOS_CATCH("");
}

void EvmKepsilonModelUtilities::CalculateTurbulentViscosityForModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY

    int number_of_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; ++i)
    {
        auto& r_node = *(rModelPart.GetCommunicator().LocalMesh().NodesBegin() + i);

        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double epsilon =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);

        const double f_mu = CalculateFmu(y_plus);

        const double nu_t = CalculateTurbulentViscosity(c_mu, tke, epsilon, f_mu);
        r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
    }

    KRATOS_CATCH("");
}

// template instantiations
template double EvmKepsilonModelUtilities::CalculateSourceTerm<2>(
    const BoundedMatrix<double, 2, 2>&, const double, const double);

template double EvmKepsilonModelUtilities::CalculateSourceTerm<3>(
    const BoundedMatrix<double, 3, 3>&, const double, const double);
///@}

} // namespace Kratos
