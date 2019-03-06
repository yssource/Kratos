#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/node.h"

#if !defined(KRATOS_RANS_APPLICATION_GEOMETRY_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_GEOMETRY_UTILITIES_H_INCLUDED

namespace Kratos
{
void CalculateGeometryData(const Geometry<Node<3>>& rGeometry,
                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                           Vector& rGaussWeights,
                           Matrix& rNContainer,
                           Geometry<Node<3>>::ShapeFunctionsGradientsType& rDN_DX);
} // namespace Kratos

#endif