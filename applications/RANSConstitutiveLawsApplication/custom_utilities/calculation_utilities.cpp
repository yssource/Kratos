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

template <class NodeType>
void LowerBound(ModelPart& rModelPart, Variable<double>& rVariable, const double MinValue)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; i++)
    {
        NodeType& r_current_node = *(rModelPart.NodesBegin() + i);
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

    // linear region
    double utau = sqrt(velocity_norm * kinematic_viscosity / wall_distance);
    double yplus = wall_distance * utau / kinematic_viscosity;

    const double limit_yplus = 11.06;
    const double inv_von_karman = 1.0 / von_karman;

    // log region
    if (yplus > limit_yplus)
    {
        // wall_vel / utau = 1/kappa * log(yplus) + B
        // this requires solving a nonlinear problem:
        // f(utau) = utau*(1/kappa * log(y*utau/nu) + B) - wall_vel = 0
        // note that f'(utau) = 1/kappa * log(y*utau/nu) + B + 1/kappa

        unsigned int iter = 0;
        double dx = 1e10;
        const double tol = 1e-6;
        double uplus = inv_von_karman * log(yplus) + beta;

        while (iter < max_iterations && fabs(dx) > tol * utau)
        {
            // Newton-Raphson iteration
            double f = utau * uplus - velocity_norm;
            double df = uplus + inv_von_karman;
            dx = f / df;

            // Update variables
            utau -= dx;
            yplus = wall_distance * utau / kinematic_viscosity;
            uplus = inv_von_karman * log(yplus) + beta;
            ++iter;
        }

        KRATOS_WARNING_IF("CalculateYplus", iter == max_iterations)
            << "Wall (logarithmic region) Newton-Raphson did not converge. "
               "residual > tolerance [ "
            << std::scientific << dx << " > " << std::scientific << tol << " ]\n";
    }

    return yplus;
}

template <class GeometryType>
double EvaluateInPoint(const GeometryType& rGeometry,
                       const Variable<double>& rVariable,
                       const Vector& rShapeFunction,
                       const int Step)
{
    const unsigned int number_of_nodes = rGeometry.PointsNumber();
    double value = 0.0;
    for (unsigned int c = 0; c < number_of_nodes; c++)
    {
        value += rShapeFunction[c] * rGeometry[c].FastGetSolutionStepValue(rVariable, Step);
    }

    return value;
}

template <class GeometryType>
array_1d<double, 3> EvaluateInPoint(const GeometryType& rGeometry,
                                    const Variable<array_1d<double, 3>>& rVariable,
                                    const Vector& rShapeFunction,
                                    const int Step)
{
    const unsigned int number_of_nodes = rGeometry.PointsNumber();
    array_1d<double, 3> value = ZeroVector(3);
    for (unsigned int c = 0; c < number_of_nodes; c++)
    {
        value += rShapeFunction[c] * rGeometry[c].FastGetSolutionStepValue(rVariable, Step);
    }

    return value;
}

template <unsigned int TDim>
double CalculateTrace(const BoundedMatrix<double, TDim, TDim>& rMatrix)
{
    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        value += rMatrix(i, i);

    return value;
}

// template instantiations
template double EvaluateInPoint<Geometry<Node<3>>>(const Geometry<Node<3>>&,
                                                   const Variable<double>&,
                                                   const Vector&,
                                                   const int);
template array_1d<double, 3> EvaluateInPoint<Geometry<Node<3>>>(
    const Geometry<Node<3>>&, const Variable<array_1d<double, 3>>&, const Vector&, const int);

template void LowerBound<Node<3>>(ModelPart&, Variable<double>&, const double);

template double CalculateTrace<2>(const BoundedMatrix<double, 2, 2>&);
template double CalculateTrace<3>(const BoundedMatrix<double, 3, 3>&);

} // namespace CalculationUtilities
} // namespace Kratos