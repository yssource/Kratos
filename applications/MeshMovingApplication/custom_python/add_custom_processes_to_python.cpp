//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//


// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/preserve_distance_process.h"


namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<PreserveDistanceProcess, PreserveDistanceProcess::Pointer, Process>(m,"PreserveDistanceProcess")
        .def(py::init<ModelPart&,ModelPart&,  Parameters>());

}

}  // namespace Python.
} // Namespace Kratos

