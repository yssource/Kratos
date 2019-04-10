//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rans_calculation_utilities.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RansVariableUtils>(m, "RansVariableUtils")
        .def(py::init<>())
        .def("AddToHistoricalNodeScalarVariable", &RansVariableUtils::AddToHistoricalNodeScalarVariable<Variable<double>>)
        ;

}

} // namespace Python.
} // Namespace Kratos
