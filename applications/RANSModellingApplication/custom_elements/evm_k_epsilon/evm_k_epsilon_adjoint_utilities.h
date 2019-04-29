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

#if !defined(KRATOS_EVM_K_EPSILON_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_K_EPSILON_ADJOINT_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "includes/ublas_interface.h"

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
} // namespace EvmKepsilonModelAdjointUtilities

///@}

} // namespace Kratos

#endif