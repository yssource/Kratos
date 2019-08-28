#if !defined(KRATOS_L2_ERROR_PROJECTION_UTILITY)
#define KRATOS_L2_ERROR_PROJECTION_UTILITY

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#include <vector>
#endif

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/qsvmsdemcoupled_data.h"
#include "includes/process_info.h"

namespace Kratos
{
class L2ErrorProjection
{

public:

    typedef Node < 3 > NodeType;
    typedef Properties PropertiesType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
    typedef ModelPart::NodesContainerType::iterator     NodeIterator;
    typedef Kratos::Vector ShapeFunctionsType;

    KRATOS_CLASS_POINTER_DEFINITION(L2ErrorProjection);

L2ErrorProjection()
{}

virtual ~L2ErrorProjection(){}

double GetL2Projection(ModelPart& r_model_part)
{
    const unsigned int n_elements = r_model_part.Elements().size();
    double area, interpolator = 0.0;
    const unsigned int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    Matrix NContainer;
    Vector N;
    std::vector<double> error;
    array_1d<double, 3> scalar_product = ZeroVector(3);

    for (unsigned int i = 0; i < n_elements; ++i){

        ElementIterator ielem = r_model_part.ElementsBegin() + i;
        GeometryType& rGeom = ielem->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();

        NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const SizeType& NumGauss = NContainer.size1();

        for (SizeType g = 0; g < NumGauss; ++g){
            N = row(NContainer, g);
        }

        area = this->GetElementArea(area, rGeom);

        for (unsigned int j = 0; j < NumNodes; ++j){
            error.push_back(rGeom[j].FastGetSolutionStepValue(ERROR_X));
            error.push_back(rGeom[j].FastGetSolutionStepValue(ERROR_Y));
            error.push_back(rGeom[j].FastGetSolutionStepValue(ERROR_Z));

            for (unsigned int k = 0; k < dim; ++k){
                for (unsigned int l = 0; l < NumNodes; ++l){
                    scalar_product[k] += error[k] * N[l];
                }
            }
            interpolator += DEM_INNER_PRODUCT_3(scalar_product, scalar_product);
            error.clear();
            scalar_product = ZeroVector(3);
        }
        interpolator *= area;
    }

    return std::sqrt(interpolator);
}

double GetElementArea(double& area, GeometryType& rGeometry)
{
    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    double x20 = rGeometry[2].X() - rGeometry[0].X();
    double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    double x30 = rGeometry[3].X() - rGeometry[0].X();
    double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

    return pow(detJ/6.0,1./3.);

}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class L2ErrorProjection

} // namespace Kratos.

#endif // KRATOS_L2_ERROR_PROJECTION_UTILITY  defined