//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <pybind11/stl.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/kutta_condition_process.h"
#include "custom_processes/define_wake_2d_process.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
	namespace py = pybind11;

    py::class_<KuttaConditionProcess, KuttaConditionProcess::Pointer, Process >
        (m, "KuttaConditionProcess")
        .def(py::init<ModelPart&>())
        ;

    py::class_<DefineWake2DProcess, DefineWake2DProcess::Pointer, Process >
        (m, "DefineWake2DProcess")
        .def(py::init<ModelPart&>());

    // py::class_<DefineWake2DProcess, DefineWake2DProcess::Pointer , Process>
    //     (m, "DefineWake2DProcess")
    //     .def(py::init<ModelPart&, Parameters>());
}

}  // namespace Python.

} // Namespace Kratos
