#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/math_utils.h"

#if !defined(KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED

namespace Kratos
{
#define CheckIfVariableIsPositive(variable) \
    KRATOS_DEBUG_ERROR_IF(variable < 0.0)   \
        << #variable << " < 0.0 [ " << std::scientific << variable << " < 0.0 ]\n";

#define PrintIfVariableIsNegative(variable)            \
    KRATOS_WARNING_IF("NegativeCheck", variable < 0.0) \
        << #variable << " < 0.0 [ " << std::scientific << variable << " < 0.0 ]\n";

#define PrintIfVariableIsPositive(variable)            \
    KRATOS_WARNING_IF("PositiveCheck", variable > 0.0) \
        << #variable << " > 0.0 [ " << std::scientific << variable << " > 0.0 ]\n";

namespace CalculationUtilities
{

bool IsPositivityPreserving(const Matrix& rA);

bool IsPositivityPreserving(const Vector& rb);

void CalculateGeometryData(const Geometry<Node<3>>& rGeometry,
                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                           Vector& rGaussWeights,
                           Matrix& rNContainer,
                           Geometry<Node<3>>::ShapeFunctionsGradientsType& rDN_DX);

Geometry<Node<3>>::ShapeFunctionsGradientsType CalculateGeometryParameterDerivatives(
    const Geometry<Node<3>>& rGeometry, const GeometryData::IntegrationMethod& rIntegrationMethod);

double CalculateYplus(const double velocity_norm,
                      const double wall_distance,
                      const double kinematic_viscosity,
                      const double von_karman,
                      const double beta,
                      const unsigned int max_iterations);

template <class NodeType>
void LowerBound(ModelPart& rModelPart, const Variable<double>& rVariable, const double MinValue);

template <class NodeType>
void UpperBound(ModelPart& rModelPart, const Variable<double>& rVariable, const double MaxValue);

template <class NodeType>
void ClipVariable(ModelPart& rModelPart,
                  const Variable<double>& rVariable,
                  const double MinValue,
                  const double MaxValue);

template <class GeometryType>
double EvaluateInPoint(const GeometryType& rGeometry,
                       const Variable<double>& rVariable,
                       const Vector& rShapeFunction,
                       const int Step = 0);

template <class GeometryType>
array_1d<double, 3> EvaluateInPoint(const GeometryType& rGeometry,
                                    const Variable<array_1d<double, 3>>& rVariable,
                                    const Vector& rShapeFunction,
                                    const int Step = 0);

template <unsigned int TDim>
double CalculateMatrixTrace(const BoundedMatrix<double, TDim, TDim>& rMatrix);

template <class NodeType>
void CheckBounds(ModelPart& rModelPart, const Variable<double>& rVariable, const std::string info = "");

template <class NodeType>
double WarnIfNegative(ModelPart& rModelPart, const Variable<double>& rVariable, const std::string info = "");

} // namespace CalculationUtilities
} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined