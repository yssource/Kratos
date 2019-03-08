#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"

#if !defined(KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED

namespace Kratos
{
#define CheckIfVariableIsPositive(variable) \
    KRATOS_DEBUG_ERROR_IF(variable < 0.0)   \
        << #variable << " < 0.0 [ " << std::scientific << variable << " < 0.0 ]\n";

namespace CalculationUtilities
{
void CalculateGeometryData(const Geometry<Node<3>>& rGeometry,
                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                           Vector& rGaussWeights,
                           Matrix& rNContainer,
                           Geometry<Node<3>>::ShapeFunctionsGradientsType& rDN_DX);

double CalculateYplus(const double velocity_norm,
                      const double wall_distance,
                      const double kinematic_viscosity,
                      const double von_karman,
                      const double beta,
                      const unsigned int max_iterations);

template <class NodeType>
void LowerBound(ModelPart& rModelPart, Variable<double>& rVariable, const double MinValue);

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
double CalculateTrace(const BoundedMatrix<double, TDim, TDim> rMatrix);

} // namespace CalculationUtilities
} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined