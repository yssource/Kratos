//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// Project includes
#include "evm_k_epsilon_adjoint_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"

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

namespace EvmKepsilonModelAdjointUtilities
{
void CalculateGaussSensitivities(Matrix& rGaussSensitivities,
                                 const Matrix& rNodalSensitivities,
                                 const Vector& rGaussShapeFunctions)
{
    const std::size_t number_of_nodes = rNodalSensitivities.size1();
    const std::size_t domain_size = rNodalSensitivities.size2();

    if (rGaussSensitivities.size1() != number_of_nodes ||
        rGaussSensitivities.size2() != domain_size)
        rGaussSensitivities.resize(number_of_nodes, domain_size);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        for (std::size_t i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            rGaussSensitivities(i_node, i_dim) =
                rGaussShapeFunctions[i_node] * rNodalSensitivities(i_node, i_dim);
        }
    }
}

void CalculateGaussSensitivities(Vector& rGaussSensitivities,
                                 const Vector& rNodalSensitivities,
                                 const Vector& rGaussShapeFunctions)
{
    const std::size_t number_of_nodes = rNodalSensitivities.size();

    if (rGaussSensitivities.size() != number_of_nodes)
        rGaussSensitivities.resize(number_of_nodes);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rGaussSensitivities[i_node] =
            rGaussShapeFunctions[i_node] * rNodalSensitivities[i_node];
    }
}

void CalculateNodalFmuVelocitySensitivities(Matrix& rFmuNodalSensitivities,
                                            const Vector& nodal_y_plus,
                                            const Matrix& rYPlusNodalSensitivities)
{
    const std::size_t number_of_nodes = rYPlusNodalSensitivities.size1();
    const std::size_t domain_size = rYPlusNodalSensitivities.size2();

    if (rFmuNodalSensitivities.size1() != number_of_nodes ||
        rFmuNodalSensitivities.size2() != domain_size)
        rFmuNodalSensitivities.resize(number_of_nodes, domain_size);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        for (std::size_t i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            rFmuNodalSensitivities(i_node, i_dim) =
                std::exp(-0.0115 * nodal_y_plus[i_node]) * 0.0115 *
                rYPlusNodalSensitivities(i_node, i_dim);
        }
    }
}

void CalculateGaussFmuVelocitySensitivities(Matrix& rFmuGaussSensitivities,
                                            const double y_plus,
                                            const Matrix& rYPlusNodalSensitivities,
                                            const Vector& rGaussShapeFunctions)
{
    CalculateGaussSensitivities(rFmuGaussSensitivities,
                                rYPlusNodalSensitivities, rGaussShapeFunctions);
    const double coeff = 0.0115 * std::exp(-0.0115 * y_plus);

    noalias(rFmuGaussSensitivities) = rFmuGaussSensitivities * coeff;
}

void CalculateNodalTurbulentViscosityVelocitySensitivities(
    Matrix& rTurbulentViscosityNodalSensitivities,
    const double c_mu,
    const Vector& nodal_turbulent_kinetic_energy,
    const Vector& nodal_turbulent_energy_dissipation_rate,
    const Matrix& rFmuNodalSensitivities)
{
    const std::size_t number_of_nodes = rFmuNodalSensitivities.size1();
    const std::size_t domain_size = rFmuNodalSensitivities.size2();

    if (rTurbulentViscosityNodalSensitivities.size1() != number_of_nodes ||
        rTurbulentViscosityNodalSensitivities.size2() != domain_size)
        rTurbulentViscosityNodalSensitivities.resize(number_of_nodes, domain_size);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        for (std::size_t i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            rTurbulentViscosityNodalSensitivities(i_node, i_dim) =
                c_mu * std::pow(nodal_turbulent_kinetic_energy[i_node], 2) /
                nodal_turbulent_energy_dissipation_rate[i_node] *
                rFmuNodalSensitivities(i_node, i_dim);
        }
    }
}

void CalculateNodalTurbulentViscosityTKESensitivities(Vector& rTurbulentViscosityNodalSensitivities,
                                                      const double c_mu,
                                                      const Vector& nodal_turbulent_kinetic_energy,
                                                      const Vector& nodal_turbulent_energy_dissipation_rate,
                                                      const Vector& nodal_f_mu)
{
    const std::size_t number_of_nodes = nodal_f_mu.size();

    if (rTurbulentViscosityNodalSensitivities.size() != number_of_nodes)
        rTurbulentViscosityNodalSensitivities.resize(number_of_nodes);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rTurbulentViscosityNodalSensitivities[i_node] =
            2 * c_mu * nodal_f_mu[i_node] * nodal_turbulent_kinetic_energy[i_node] /
            nodal_turbulent_energy_dissipation_rate[i_node];
    }
}

void CalculateNodalTurbulentViscosityEpsilonSensitivities(
    Vector& rTurbulentViscosityNodalSensitivities,
    const double c_mu,
    const Vector& nodal_turbulent_kinetic_energy,
    const Vector& nodal_turbulent_energy_dissipation_rate,
    const Vector& nodal_f_mu)
{
    const std::size_t number_of_nodes = nodal_f_mu.size();

    if (rTurbulentViscosityNodalSensitivities.size() != number_of_nodes)
        rTurbulentViscosityNodalSensitivities.resize(number_of_nodes);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rTurbulentViscosityNodalSensitivities[i_node] =
            -1.0 * c_mu * nodal_f_mu[i_node] *
            std::pow(nodal_turbulent_kinetic_energy[i_node] /
                         nodal_turbulent_energy_dissipation_rate[i_node],
                     2);
    }
}
template <unsigned int TDim>
void CalculateProductionVelocitySensitivities(Matrix& rOutput,
                                              const double NuT,
                                              const Matrix& rNuTVelocityDerivatives,
                                              const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                                              const Matrix& rShapeDerivatives)
{
    const std::size_t number_of_nodes = rShapeDerivatives.size1();

    if (rOutput.size1() != number_of_nodes || rOutput.size2() != TDim)
        rOutput.resize(number_of_nodes, TDim);
    rOutput.clear();

    double velocity_divergence = 0.0;
    velocity_divergence =
        RansCalculationUtilities().CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) = rVelocityGradient + trans(rVelocityGradient) -
                                      (2.0 / 3.0) * velocity_divergence * identity;

    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            value += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    noalias(rOutput) = rNuTVelocityDerivatives * value;

    for (std::size_t c = 0; c < number_of_nodes; ++c)
    {
        for (std::size_t k = 0; k < TDim; ++k)
        {
            value = 0.0;

            for (std::size_t j = 0; j < TDim; ++j)
            {
                value += rShapeDerivatives(c, j) * rVelocityGradient(k, j);
                value += rShapeDerivatives(c, j) * rVelocityGradient(j, k);
                value -= rShapeDerivatives(c, k) * (2.0 / 3.0) * rVelocityGradient(j, j);
                value += rShapeDerivatives(c, j) * reynolds_stress_tensor(k, j);
            }

            rOutput(c, k) += NuT * value;
        }
    }
}

template <unsigned int TDim>
void CalculateProductionScalarSensitivities(Vector& rOutput,
                                            const Vector& rNuTScalarDerivatives,
                                            const BoundedMatrix<double, TDim, TDim>& rVelocityGradient)
{
    const std::size_t number_of_nodes = rNuTScalarDerivatives.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);

    double velocity_divergence = 0.0;
    velocity_divergence =
        RansCalculationUtilities().CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) = rVelocityGradient + trans(rVelocityGradient) -
                                      (2.0 / 3.0) * velocity_divergence * identity;

    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            value += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    noalias(rOutput) = rNuTScalarDerivatives * value;
}

void CalculateThetaVelocitySensitivity(Matrix& rOutput,
                                       const double c_mu,
                                       const double f_mu,
                                       const double tke,
                                       const double nu_t,
                                       const Matrix& rFmuSensitivities,
                                       const Matrix& rNuTSensitivities)
{
    std::size_t number_of_nodes = rFmuSensitivities.size1();
    std::size_t domain_size = rFmuSensitivities.size2();

    if (rOutput.size1() != number_of_nodes || rOutput.size2() != domain_size)
        rOutput.resize(number_of_nodes, domain_size);
    rOutput.clear();

    noalias(rOutput) += rFmuSensitivities * c_mu * tke / nu_t;
    noalias(rOutput) -= rNuTSensitivities * c_mu * f_mu * tke / std::pow(nu_t, 2);
}

void CalculateThetaTKESensitivity(Vector& rOutput,
                                  const double c_mu,
                                  const double f_mu,
                                  const double tke,
                                  const double nu_t,
                                  const Vector& rNuTGaussSensitivities,
                                  const Vector& rGaussShapeFunctions)
{
    std::size_t number_of_nodes = rNuTGaussSensitivities.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);
    rOutput.clear();

    noalias(rOutput) += rGaussShapeFunctions * c_mu * f_mu / nu_t;
    noalias(rOutput) -= rNuTGaussSensitivities * c_mu * f_mu * tke / std::pow(nu_t, 2);
}

void CalculateThetaEpsilonSensitivity(Vector& rOutput,
                                      const double c_mu,
                                      const double f_mu,
                                      const double tke,
                                      const double nu_t,
                                      const Vector& rNuTSensitivities)
{
    std::size_t number_of_nodes = rNuTSensitivities.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);

    noalias(rOutput) =
        rNuTSensitivities * (-1.0 * c_mu * f_mu * tke / std::pow(nu_t, 2));
}

void CalculateTurbulentReynoldsNumberVelocitySensitivity(Matrix& rOutput,
                                                         const double tke,
                                                         const double epsilon,
                                                         const double nu_t,
                                                         const Matrix& rNuTSensitivities)
{
    std::size_t number_of_nodes = rNuTSensitivities.size1();
    std::size_t domain_size = rNuTSensitivities.size2();

    if (rOutput.size1() != number_of_nodes || rOutput.size2() != domain_size)
        rOutput.resize(number_of_nodes, domain_size);

    noalias(rOutput) = rNuTSensitivities * (-1.0 * std::pow(tke / nu_t, 2) / epsilon);
}

void CalculateTurbulentReynoldsNumberTKESensitivity(Vector& rOutput,
                                                    const double tke,
                                                    const double epsilon,
                                                    const double nu_t,
                                                    const Vector& rNuTSensitivities,
                                                    const Vector& rGaussShapeFunctions)
{
    std::size_t number_of_nodes = rNuTSensitivities.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);

    noalias(rOutput) = rGaussShapeFunctions * (2.0 * tke / (epsilon * nu_t));
    noalias(rOutput) -= rNuTSensitivities * (std::pow(tke / nu_t, 2) / epsilon);
}

void CalculateTurbulentReynoldsNumberEpsilonSensitivity(Vector& rOutput,
                                                        const double tke,
                                                        const double epsilon,
                                                        const double nu_t,
                                                        const Vector& rNuTSensitivities,
                                                        const Vector& rGaussShapeFunctions)
{
    std::size_t number_of_nodes = rNuTSensitivities.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);

    noalias(rOutput) = rGaussShapeFunctions * (-1.0 * std::pow(tke / epsilon, 2) / nu_t);
    noalias(rOutput) -= rNuTSensitivities * (std::pow(tke / nu_t, 2) / epsilon);
}

void CalculateF2VelocitySensitivity(Matrix& rOutput,
                                    const double tke,
                                    const double epsilon,
                                    const double nu_t,
                                    const Matrix& rNuTSensitivities)
{
    CalculateTurbulentReynoldsNumberVelocitySensitivity(
        rOutput, tke, epsilon, nu_t, rNuTSensitivities);

    if (std::abs(epsilon) <= std::numeric_limits<double>::epsilon())
    {
        rOutput.clear();
        return;
    }

    const double Re_t = std::pow(tke, 2) / (nu_t * epsilon * 6.0);
    const double coeff = 0.44 * Re_t * std::exp(-std::pow(Re_t, 2));

    noalias(rOutput) = rOutput * coeff;
}

void CalculateF2ScalarSensitivity(Vector& rOutput,
                                  const double epsilon,
                                  const double Re_t,
                                  const Vector& rReTSensitivities)
{
    std::size_t number_of_nodes = rReTSensitivities.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);

    if (std::abs(epsilon) <= std::numeric_limits<double>::epsilon())
    {
        rOutput.clear();
        return;
    }
    noalias(rOutput) = rReTSensitivities *
                       (0.44 * Re_t / 6.0 * std::exp(-std::pow(Re_t / 6.0, 2)));
}

void CalculateReactionVelocitySensitivities(Matrix& rOutput,
                                            const double wall_distance,
                                            const Matrix& rNuTSensitivities,
                                            const Matrix& rThetaSensitivities)
{
    std::size_t number_of_nodes = rNuTSensitivities.size1();
    std::size_t domain_size = rNuTSensitivities.size2();

    if (rOutput.size1() != number_of_nodes || rOutput.size2() != domain_size)
        rOutput.resize(number_of_nodes, domain_size);

    noalias(rOutput) = rNuTSensitivities * (2.0 / std::pow(wall_distance, 2)) + rThetaSensitivities;
}

void CalculateReactionScalarSensitivities(Vector& rOutput,
                                          const double wall_distance,
                                          const Vector& rNuTSensitivities,
                                          const Vector& rThetaSensitivities)
{
    std::size_t number_of_nodes = rNuTSensitivities.size();

    if (rOutput.size() != number_of_nodes)
        rOutput.resize(number_of_nodes);

    noalias(rOutput) = rNuTSensitivities * (2.0 / std::pow(wall_distance, 2)) + rThetaSensitivities;
}

template void CalculateProductionVelocitySensitivities<2>(
    Matrix&, const double, const Matrix&, const BoundedMatrix<double, 2, 2>&, const Matrix&);
template void CalculateProductionVelocitySensitivities<3>(
    Matrix&, const double, const Matrix&, const BoundedMatrix<double, 3, 3>&, const Matrix&);

template void CalculateProductionScalarSensitivities<2>(
    Vector&, const Vector&, const BoundedMatrix<double, 2, 2>&);

template void CalculateProductionScalarSensitivities<3>(
    Vector&, const Vector&, const BoundedMatrix<double, 3, 3>&);

} // namespace EvmKepsilonModelAdjointUtilities

///@}

} // namespace Kratos
