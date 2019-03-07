#include "calculation_utilities.h"
#include <cmath>

namespace Kratos
{
namespace CalculationUtilities
{
void CalculateGeometryData(const Geometry<Node<3>>& rGeometry,
                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                           Vector& rGaussWeights,
                           Matrix& rNContainer,
                           Geometry<Node<3>>::ShapeFunctionsGradientsType& rDN_DX)
{
    const unsigned int number_of_gauss_points =
        rGeometry.IntegrationPointsNumber(rIntegrationMethod);

    Vector DetJ;
    rGeometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, rIntegrationMethod);

    const std::size_t number_of_nodes = rGeometry.PointsNumber();

    if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != number_of_nodes)
    {
        rNContainer.resize(number_of_gauss_points, number_of_nodes, false);
    }
    rNContainer = rGeometry.ShapeFunctionsValues(rIntegrationMethod);

    const Geometry<Node<3>>::IntegrationPointsArrayType& IntegrationPoints =
        rGeometry.IntegrationPoints(rIntegrationMethod);

    if (rGaussWeights.size() != number_of_gauss_points)
    {
        rGaussWeights.resize(number_of_gauss_points, false);
    }

    for (unsigned int g = 0; g < number_of_gauss_points; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

void LowerBound(ModelPart& rModelPart, Variable<double>& rVariable, const double MinValue)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; i++)
    {
        Node<3>& r_current_node = *(rModelPart.NodesBegin() + i);
        double& value = r_current_node.FastGetSolutionStepValue(rVariable);
        value = std::max<double>(MinValue, value);
    }
}

double CalculateYplus(const double velocity_norm,
                      const double wall_distance,
                      const double kinematic_viscosity,
                      const double von_karman,
                      const double beta,
                      const unsigned int max_iterations)
{
    CheckIfVariableIsPositive(velocity_norm);
    CheckIfVariableIsPositive(wall_distance);
    CheckIfVariableIsPositive(kinematic_viscosity);
    CheckIfVariableIsPositive(von_karman);
    CheckIfVariableIsPositive(beta);

    // try linear law
    double y_plus = std::sqrt(velocity_norm * wall_distance / kinematic_viscosity);

    // If the linear low doesnt match within the range, try logrithmic law
    if (y_plus > 11.06)
    {
        unsigned int i;
        double u_tau = std::sqrt(velocity_norm * kinematic_viscosity / wall_distance);
        double prev_u_tau = 0.0;
        for (i = 0; i < max_iterations; ++i)
        {
            prev_u_tau = u_tau;
            u_tau = velocity_norm /
                    (std::log(u_tau * wall_distance / kinematic_viscosity) / von_karman + beta);
        }
        const double delta_u_tau = std::abs(u_tau - prev_u_tau);
        KRATOS_INFO_IF("TurbulenceEvmProcess", delta_u_tau > 1e-5)
            << "WARNING: Maximum number of iterations reached for y_plus "
               "calculation. error_u_tau = "
            << std::scientific << delta_u_tau << ".\n";
        y_plus = u_tau * wall_distance / kinematic_viscosity;
    }
    return y_plus;
}
} // namespace CalculationUtilities
} // namespace Kratos