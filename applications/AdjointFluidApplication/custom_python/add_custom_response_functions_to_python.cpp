// External includes
#include <boost/python.hpp>

// Project includes

// Application includes
#include "custom_utilities/response_function.h"
#include "custom_utilities/drag_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"
#include "custom_utilities/potentialflow_lift_response_function.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void CalculateSensitivityGradient1(ResponseFunction& dummy,Element& rAdjointElem,
  const Variable<array_1d<double,3>>& rVariable,
  const Matrix& rDerivativesMatrix,
  Vector& rResponseGradient,
  ProcessInfo& rProcessInfo)
{
dummy.CalculateSensitivityGradient(rAdjointElem,rVariable,rDerivativesMatrix,rResponseGradient,rProcessInfo);
}

void CalculateSensitivityGradient2(ResponseFunction& dummy,Element& rAdjointElem,
  const Variable<double>& rVariable,
  const Matrix& rDerivativesMatrix,
  Vector& rResponseGradient,
  ProcessInfo& rProcessInfo)
{
dummy.CalculateSensitivityGradient(rAdjointElem,rVariable,rDerivativesMatrix,rResponseGradient,rProcessInfo);
}

void AddCustomResponseFunctionsToPython()
{
  class_<ResponseFunction, boost::noncopyable>("ResponseFunction", init<ModelPart&, Parameters&>())
        .def("Initialize", &ResponseFunction::Initialize)
        .def("InitializeSolutionStep", &ResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ResponseFunction::FinalizeSolutionStep)
        .def("Check", &ResponseFunction::Check)
        .def("Clear", &ResponseFunction::Clear)
        .def("CalculateGradient",
             &ResponseFunction::CalculateGradient)
        .def("CalculateFirstDerivativesGradient",
             &ResponseFunction::CalculateFirstDerivativesGradient)
        .def("CalculateSecondDerivativesGradient",
             &ResponseFunction::CalculateSecondDerivativesGradient)
        .def("CalculateSensitivityGradient", &CalculateSensitivityGradient1)
        .def("CalculateSensitivityGradient", &CalculateSensitivityGradient2)         
        .def("CalculateValue", &ResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &ResponseFunction::UpdateSensitivities);

    class_<DragResponseFunction<2>, bases<ResponseFunction>, boost::noncopyable>
      ("DragResponseFunction2D", init<ModelPart&, Parameters&>());

    class_<DragResponseFunction<3>, bases<ResponseFunction>, boost::noncopyable>
      ("DragResponseFunction3D", init<ModelPart&, Parameters&>());

    class_<PotentialFlowLiftResponseFunction<2>, bases<ResponseFunction>, boost::noncopyable>
      ("PotentialFlowLiftResponseFunction2D", init<ModelPart&, Parameters&>());

    class_<PotentialFlowLiftResponseFunction<3>, bases<ResponseFunction>, boost::noncopyable>
      ("PotentialFlowLiftResponseFunction3D", init<ModelPart&, Parameters&>());

}

} // namespace Python

} // namespace Kratos
