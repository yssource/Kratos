#include "rans_calculation_utilities.h"

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

void RansCalculationUtilities::CalculateGeometryParameterDerivativesShapeSensitivity(
    Matrix& rOutput, const ShapeParameter& rShapeDerivative, const Matrix& rDnDe, const Matrix& rDeDx)
{
    std::size_t domain_size = rDeDx.size1();
    if (rOutput.size1() != domain_size || rOutput.size2() != domain_size)
        rOutput.resize(domain_size, domain_size);

    const Vector& r_dnc_de = row(rDnDe, rShapeDerivative.NodeIndex);

    for (std::size_t j = 0; j < domain_size; ++j)
    {
        const Vector& r_de_dxj = column(rDeDx, j);
        for (std::size_t i = 0; i < domain_size; ++i)
        {
            rOutput(i, j) = -1.0 * rDeDx(i, rShapeDerivative.Direction) *
                            inner_prod(r_dnc_de, r_de_dxj);
        }
    }
}

template <unsigned int TDim>
void RansCalculationUtilities::CalculateGradient(BoundedMatrix<double, TDim, TDim>& rOutput,
                                                 const Geometry<ModelPart::NodeType>& rGeometry,
                                                 const Variable<array_1d<double, 3>>& rVariable,
                                                 const Matrix& rShapeDerivatives,
                                                 const int Step) const
{
    rOutput.clear();
    std::size_t number_of_nodes = rGeometry.PointsNumber();

    for (unsigned int a = 0; a < number_of_nodes; ++a)
    {
        const array_1d<double, 3>& r_value =
            rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (unsigned int i = 0; i < TDim; ++i)
        {
            for (unsigned int j = 0; j < TDim; ++j)
            {
                rOutput(i, j) += rShapeDerivatives(a, j) * r_value[i];
            }
        }
    }
}

void RansCalculationUtilities::CalculateGradient(array_1d<double, 3>& rOutput,
                                                 const Geometry<ModelPart::NodeType>& rGeometry,
                                                 const Variable<double>& rVariable,
                                                 const Matrix& rShapeDerivatives,
                                                 const int Step) const
{
    rOutput.clear();
    std::size_t number_of_nodes = rGeometry.PointsNumber();
    unsigned int domain_size = rShapeDerivatives.size2();

    for (std::size_t a = 0; a < number_of_nodes; ++a)
    {
        const double value = rGeometry[a].FastGetSolutionStepValue(rVariable, Step);
        for (unsigned int i = 0; i < domain_size; ++i)
            rOutput[i] += rShapeDerivatives(a, i) * value;
    }
}

template <unsigned int TDim>
void RansCalculationUtilities::CalculateVelocityGradientSensitivities(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const int VelocityDerivNodeIndex,
    const int VelocityDerivDirection,
    const Matrix& rShapeDerivatives) const
{
    rOutput.clear();
    for (unsigned int j = 0; j < TDim; ++j)
    {
        rOutput(VelocityDerivDirection, j) =
            rShapeDerivatives(VelocityDerivNodeIndex, j);
    }
}

template <unsigned int TDim>
Vector RansCalculationUtilities::GetVector(const array_1d<double, 3>& rVector) const
{
    Vector result(TDim);

    for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
        result[i_dim] = rVector[i_dim];

    return result;
}

Vector RansCalculationUtilities::GetVector(const array_1d<double, 3>& rVector,
                                           const unsigned int Dim) const
{
    Vector result(Dim);

    for (unsigned int i_dim = 0; i_dim < Dim; ++i_dim)
        result[i_dim] = rVector[i_dim];

    return result;
}

void RansCalculationUtilities::AssembleElementVector(Vector& rElementVector,
                                                     const Vector& rSubVector,
                                                     const std::size_t NumberOfNodes,
                                                     const std::size_t Offset)
{
    const std::size_t element_local_size = rElementVector.size() / NumberOfNodes;
    const std::size_t sub_vector_local_size = rSubVector.size() / NumberOfNodes;

    for (std::size_t i = 0; i < rSubVector.size(); ++i)
    {
        const std::size_t i_node = i / sub_vector_local_size;
        const std::size_t i_dim = i % sub_vector_local_size;
        rElementVector[i_node * element_local_size + Offset + i_dim] = rSubVector[i];
    }
}

void RansCalculationUtilities::AssembleElementMatrix(Matrix& rElementMatrix,
                                                     const Matrix& rSubMatrix,
                                                     const std::size_t NumberOfNodes,
                                                     const std::size_t DerivativeOffset,
                                                     const std::size_t EquationOffset)
{
    const std::size_t element_matrix_eq_local_size = rElementMatrix.size2() / NumberOfNodes;
    const std::size_t element_matrix_deriv_local_size = rElementMatrix.size1() / NumberOfNodes;

    const std::size_t sub_matrix_derivatives_size = rSubMatrix.size1();
    const std::size_t sub_matrix_equations_size = rSubMatrix.size2();
    const std::size_t sub_matrix_equation_local_size = sub_matrix_equations_size / NumberOfNodes;
    const std::size_t sub_matrix_derivatives_local_size =
        sub_matrix_derivatives_size / NumberOfNodes;

    for (std::size_t i_deriv = 0; i_deriv < sub_matrix_derivatives_size; ++i_deriv)
    {
        const std::size_t node_deriv = i_deriv / sub_matrix_derivatives_local_size;
        const std::size_t node_deriv_dim = i_deriv % sub_matrix_derivatives_local_size;

        for (std::size_t i_eq = 0; i_eq < sub_matrix_equations_size; ++i_eq)
        {
            const std::size_t node_eq = i_eq / sub_matrix_equation_local_size;
            const std::size_t node_eq_dim = i_eq % sub_matrix_equation_local_size;

            const std::size_t element_eq =
                node_eq * element_matrix_eq_local_size + EquationOffset + node_eq_dim;
            const std::size_t element_deriv = node_deriv * element_matrix_deriv_local_size +
                                              DerivativeOffset + node_deriv_dim;

            rElementMatrix(element_deriv, element_eq) = rSubMatrix(i_deriv, i_eq);
        }
    }
}

void RansCalculationUtilities::GetSubMatrix(Matrix& rSubMatrix,
                                            const Matrix& rElementMatrix,
                                            const std::size_t NumberOfNodes,
                                            const std::size_t NodeDerivativeOffset,
                                            const std::size_t NodeEquationOffset,
                                            const std::size_t NumberOfNodalDerivatives,
                                            const std::size_t NumberOfNodalEquations)
{
    const std::size_t sub_matrix_equations = NumberOfNodes * NumberOfNodalEquations;
    const std::size_t sub_matrix_derivatives = NumberOfNodes * NumberOfNodalDerivatives;
    const std::size_t element_matrix_deriv_local_size = rElementMatrix.size1() / NumberOfNodes;
    const std::size_t element_matrix_eq_local_size = rElementMatrix.size2() / NumberOfNodes;

    if (rSubMatrix.size1() != sub_matrix_derivatives || rSubMatrix.size2() != sub_matrix_equations)
        rSubMatrix.resize(sub_matrix_derivatives, sub_matrix_equations, false);

    for (std::size_t i_deriv_node = 0; i_deriv_node < NumberOfNodes; ++i_deriv_node)
    {
        for (std::size_t i_deriv_dim = 0; i_deriv_dim < NumberOfNodalDerivatives; ++i_deriv_dim)
        {
            for (std::size_t i_eq_node = 0; i_eq_node < NumberOfNodes; ++i_eq_node)
            {
                for (std::size_t i_eq_dim = 0; i_eq_dim < NumberOfNodalEquations; ++i_eq_dim)
                {
                    const std::size_t s_i =
                        i_deriv_node * NumberOfNodalDerivatives + i_deriv_dim;
                    const std::size_t s_j = i_eq_node * NumberOfNodalEquations + i_eq_dim;
                    const std::size_t g_i = i_deriv_node * element_matrix_deriv_local_size +
                                            NodeDerivativeOffset + i_deriv_dim;
                    const std::size_t g_j = i_eq_node * element_matrix_eq_local_size +
                                            NodeEquationOffset + i_eq_dim;
                    rSubMatrix(s_i, s_j) = rElementMatrix(g_i, g_j);
                }
            }
        }
    }
}

// template instantiations

template double RansCalculationUtilities::CalculateMatrixTrace<2>(
    const BoundedMatrix<double, 2, 2>&);
template double RansCalculationUtilities::CalculateMatrixTrace<3>(
    const BoundedMatrix<double, 3, 3>&);

template void RansCalculationUtilities::CalculateGradient<2>(
    BoundedMatrix<double, 2, 2>&,
    const Geometry<ModelPart::NodeType>&,
    const Variable<array_1d<double, 3>>&,
    const Matrix&,
    const int) const;

template void RansCalculationUtilities::CalculateGradient<3>(
    BoundedMatrix<double, 3, 3>&,
    const Geometry<ModelPart::NodeType>&,
    const Variable<array_1d<double, 3>>&,
    const Matrix&,
    const int) const;

template void RansCalculationUtilities::CalculateVelocityGradientSensitivities<2>(
    BoundedMatrix<double, 2, 2>&, const int, const int, const Matrix&) const;

template void RansCalculationUtilities::CalculateVelocityGradientSensitivities<3>(
    BoundedMatrix<double, 3, 3>&, const int, const int, const Matrix&) const;

template Vector RansCalculationUtilities::GetVector<2>(const array_1d<double, 3>&) const;
template Vector RansCalculationUtilities::GetVector<3>(const array_1d<double, 3>&) const;

///@}

} // namespace Kratos