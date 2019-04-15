#include "rans_calculation_utilities.h"
#include <cmath>

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

void RansCalculationUtilities::CalculateGeometryData(
    const RansCalculationUtilities::GeometryType& rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod,
    Vector& rGaussWeights,
    Matrix& rNContainer,
    RansCalculationUtilities::GeometryType::ShapeFunctionsGradientsType& rDN_DX)
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

    const RansCalculationUtilities::GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeometry.IntegrationPoints(rIntegrationMethod);

    if (rGaussWeights.size() != number_of_gauss_points)
    {
        rGaussWeights.resize(number_of_gauss_points, false);
    }

    for (unsigned int g = 0; g < number_of_gauss_points; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

void RansCalculationUtilities::LowerBound(ModelPart& rModelPart,
                                          const Variable<double>& rVariable,
                                          const double MinValue)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    unsigned int count_of_nodes = 0;

#pragma omp parallel for reduction(+ : count_of_nodes)
    for (int i = 0; i < number_of_nodes; i++)
    {
        NodeType& r_current_node = *(rModelPart.NodesBegin() + i);
        double& value = r_current_node.FastGetSolutionStepValue(rVariable);
        if (value < MinValue)
        {
            value = MinValue;
            count_of_nodes++;
        }
    }

    KRATOS_WARNING_IF("LowerBound", count_of_nodes > 0 && mEchoLevel > 0)
        << rVariable.Name() << " of " << count_of_nodes << " nodes are less than "
        << std::scientific << MinValue << " out of total number of "
        << number_of_nodes << " nodes in " << rModelPart.Name() << ".\n";
}

void RansCalculationUtilities::ClipVariable(ModelPart& rModelPart,
                                            const Variable<double>& rVariable,
                                            const double MinValue,
                                            const double MaxValue)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    unsigned int count_of_nodes_negative = 0;
    unsigned int count_of_nodes_positive = 0;

#pragma omp parallel for reduction(+ : count_of_nodes_negative, count_of_nodes_positive)
    for (int i = 0; i < number_of_nodes; i++)
    {
        NodeType& r_current_node = *(rModelPart.NodesBegin() + i);
        double& value = r_current_node.FastGetSolutionStepValue(rVariable);
        if (value < MinValue)
        {
            value = MinValue;
            count_of_nodes_negative++;
        }
        else if (value > MaxValue)
        {
            value = MaxValue;
            count_of_nodes_positive++;
        }
    }

    const unsigned int count_of_nodes = count_of_nodes_negative + count_of_nodes_positive;

    KRATOS_WARNING_IF("ClipVariable", count_of_nodes > 0)
        << rVariable.Name() << " of " << count_of_nodes << " nodes are not in the range of [ "
        << std::scientific << MinValue << ", " << std::scientific << MaxValue
        << " ] out of total number of " << number_of_nodes << " nodes in "
        << rModelPart.Name() << " [ " << count_of_nodes_negative << " nodes < "
        << std::scientific << MinValue << "; " << count_of_nodes_positive
        << " nodes > " << std::scientific << MaxValue << " ].\n";
}

void RansCalculationUtilities::UpperBound(ModelPart& rModelPart,
                                          const Variable<double>& rVariable,
                                          const double MaxValue)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    unsigned int count_of_nodes = 0;

#pragma omp parallel for reduction(+ : count_of_nodes)
    for (int i = 0; i < number_of_nodes; i++)
    {
        NodeType& r_current_node = *(rModelPart.NodesBegin() + i);
        double& value = r_current_node.FastGetSolutionStepValue(rVariable);
        if (value > MaxValue)
        {
            value = MaxValue;
            count_of_nodes++;
        }
    }

    KRATOS_WARNING_IF("UpperBound", count_of_nodes > 0 && mEchoLevel > 0)
        << rVariable.Name() << " of " << count_of_nodes << " nodes are greater than "
        << std::scientific << MaxValue << " out of total number of "
        << number_of_nodes << " nodes in " << rModelPart.Name() << ".\n";
}

double RansCalculationUtilities::CalculateYplus(const double velocity_norm,
                                                const double wall_distance,
                                                const double kinematic_viscosity,
                                                const double von_karman,
                                                const double beta,
                                                const unsigned int max_iterations)
{
    // linear region
    double utau = sqrt(velocity_norm * kinematic_viscosity / wall_distance);
    double yplus = wall_distance * utau / kinematic_viscosity;

    const double limit_yplus = 11.06;
    const double inv_von_karman = 1.0 / von_karman;

    // log region
    if (yplus > limit_yplus)
    {
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

    KRATOS_WARNING_IF("CalculateYplus", yplus < 0.0)
        << "Calculated y_plus < 0.0 [ " << std::scientific << yplus << " < 0.0 ]\n";

    return yplus;
}

double RansCalculationUtilities::EvaluateInPoint(const RansCalculationUtilities::GeometryType& rGeometry,
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

array_1d<double, 3> RansCalculationUtilities::EvaluateInPoint(
    const RansCalculationUtilities::GeometryType& rGeometry,
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
double RansCalculationUtilities::CalculateMatrixTrace(const BoundedMatrix<double, TDim, TDim>& rMatrix)
{
    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        value += rMatrix(i, i);

    return value;
}

RansCalculationUtilities::GeometryType::ShapeFunctionsGradientsType RansCalculationUtilities::CalculateGeometryParameterDerivatives(
    const RansCalculationUtilities::GeometryType& rGeometry,
    const GeometryData::IntegrationMethod& rIntegrationMethod)
{
    const RansCalculationUtilities::GeometryType::ShapeFunctionsGradientsType& DN_De =
        rGeometry.ShapeFunctionsLocalGradients(rIntegrationMethod);
    const std::size_t number_of_nodes = rGeometry.PointsNumber();
    const unsigned int number_of_gauss_points =
        rGeometry.IntegrationPointsNumber(rIntegrationMethod);
    const std::size_t dim = rGeometry.WorkingSpaceDimension();

    RansCalculationUtilities::GeometryType::ShapeFunctionsGradientsType de_dx(number_of_gauss_points);

    Matrix geometry_coordinates(dim, number_of_nodes);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const array_1d<double, 3>& r_coordinates =
            rGeometry.Points()[i_node].Coordinates();
        for (std::size_t d = 0; d < dim; ++d)
            geometry_coordinates(d, i_node) = r_coordinates[d];
    }

    for (unsigned int g = 0; g < number_of_gauss_points; ++g)
    {
        const Matrix& r_current_local_gradients = DN_De[g];
        Matrix current_dx_de(dim, dim);
        noalias(current_dx_de) = prod(geometry_coordinates, r_current_local_gradients);
        Matrix inv_current_dx_de(dim, dim);
        double det_J;
        MathUtils<double>::InvertMatrix<Matrix, Matrix>(
            current_dx_de, inv_current_dx_de, det_J);

        de_dx[g] = inv_current_dx_de;
    }

    return de_dx;
}

void RansCalculationUtilities::CheckBounds(ModelPart& rModelPart,
                                           const Variable<double>& rVariable,
                                           const std::string info)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    double min{0.0}, max{0.0};
    bool initialized = false;

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; i++)
    {
        NodeType& r_current_node = *(rModelPart.NodesBegin() + i);
        const double value = r_current_node.FastGetSolutionStepValue(rVariable);

        if (!initialized)
        {
            min = value;
            max = value;
            initialized = true;
        }

        min = std::min(min, value);
        max = std::max(max, value);
    }

    KRATOS_INFO("CheckBounds")
        << info << " - " << rVariable.Name() << " is bounded [ "
        << std::scientific << min << ", " << std::scientific << max << " ].\n";
}

double RansCalculationUtilities::WarnIfNegative(ModelPart& rModelPart,
                                                const Variable<double>& rVariable,
                                                const std::string info)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    unsigned int count_of_nodes = 0;

    double aggregated_negatives = 0.0;

#pragma omp parallel for reduction(+ : count_of_nodes, aggregated_negatives)
    for (int i = 0; i < number_of_nodes; i++)
    {
        NodeType& r_current_node = *(rModelPart.NodesBegin() + i);
        const double value = r_current_node.FastGetSolutionStepValue(rVariable);
        if (value < 0.0)
        {
            count_of_nodes++;
            aggregated_negatives += value;
        }
    }

    KRATOS_WARNING_IF("WarnIfNegative", count_of_nodes > 0)
        << rVariable.Name() << " of " << count_of_nodes << " nodes are less than "
        << std::scientific << 0.0 << " out of total number of "
        << number_of_nodes << " nodes in " << rModelPart.Name() << ".\n";

    return aggregated_negatives;
}

// template instantiations

template double RansCalculationUtilities::CalculateMatrixTrace<2>(
    const BoundedMatrix<double, 2, 2>&);
template double RansCalculationUtilities::CalculateMatrixTrace<3>(
    const BoundedMatrix<double, 3, 3>&);

///@}

} // namespace Kratos