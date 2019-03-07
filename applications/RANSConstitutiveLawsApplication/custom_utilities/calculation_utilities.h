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

void LowerBound(ModelPart& rModelPart, Variable<double>& rVariable, const double MinValue);

double CalculateYplus(const double velocity_norm,
                      const double wall_distance,
                      const double kinematic_viscosity,
                      const double von_karman,
                      const double beta,
                      const unsigned int max_iterations);

} // namespace CalculationUtilities
} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined