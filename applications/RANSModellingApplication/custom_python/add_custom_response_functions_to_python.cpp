// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_response_functions/rans_drag_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{

template<unsigned TDim>
class PyRansDragResponseFunction : public RansDragResponseFunction<TDim>
{
    public:
        using RansDragResponseFunction<TDim>::RansDragResponseFunction;
        void Initialize() override {
            PYBIND11_OVERLOAD(void, RansDragResponseFunction<TDim>, Initialize, );
        }
        void InitializeSolutionStep() override {
            PYBIND11_OVERLOAD(void, RansDragResponseFunction<TDim>, InitializeSolutionStep, );
        }
        void FinalizeSolutionStep() override {
            PYBIND11_OVERLOAD(void, RansDragResponseFunction<TDim>, FinalizeSolutionStep, );
        }
        double CalculateValue(ModelPart& rModelPart) override {
            PYBIND11_OVERLOAD_PURE(double, RansDragResponseFunction<TDim>, CalculateValue, rModelPart);
        }
};


namespace Python
{

void AddCustomResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;
    py::class_<
        RansDragResponseFunction<2>,
        PyRansDragResponseFunction<2>,
        RansDragResponseFunction<2>::Pointer,
        AdjointResponseFunction>(m,"RansDragResponseFunction2D")
        .def(py::init<Parameters, ModelPart&>());

    py::class_<
        RansDragResponseFunction<3>,
        PyRansDragResponseFunction<3>,
        RansDragResponseFunction<3>::Pointer,
        AdjointResponseFunction>(m,"RansDragResponseFunction3D")
        .def(py::init<Parameters, ModelPart&>());

}

} // namespace Python

} // namespace Kratos
