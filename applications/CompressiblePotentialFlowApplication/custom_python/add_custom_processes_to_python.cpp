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


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/kutta_condition_process.h"
#include "custom_processes/compute_lift_level_set_process.h"
#include "custom_processes/metrics_potential_levelset_process.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
	namespace py = pybind11;

    py::class_<KuttaConditionProcess, KuttaConditionProcess::Pointer, Process >
        (m, "KuttaConditionProcess")
        .def(py::init<ModelPart&>())
        ;

    py::class_<ComputeLiftLevelSetProcess, ComputeLiftLevelSetProcess::Pointer, Process >
        (m, "ComputeLiftLevelSetProcess")
        .def(py::init<ModelPart&,Vector&>())
        ;

    py::class_<ComputePotentialLevelSetSolMetricProcess<2>, ComputePotentialLevelSetSolMetricProcess<2>::Pointer, Process>
        (m, "ComputePotentialLevelSetSolMetricProcess2D")
        .def(py::init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(py::init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        ;
}

}
} // Namespace Kratos
