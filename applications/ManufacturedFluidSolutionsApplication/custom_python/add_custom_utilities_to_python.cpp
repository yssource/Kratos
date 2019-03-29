//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes

// External includes

// Project includes
#include "add_custom_utilities_to_python.h"
#include "custom_utilities/manufactured_solution_utility.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

     py::class_< ManufacturedSolutionUtility, ManufacturedSolutionUtility::Pointer>(m, "ManufacturedSolutionUtility")
     	.def(py::init<ModelPart&, ManufacturedSolution&>() )
     	.def("SetBodyForce", &ManufacturedSolutionUtility::SetBodyForce)
     	.def("SetVelocity", &ManufacturedSolutionUtility::SetVelocity)
     	.def("SetPressure", &ManufacturedSolutionUtility::SetPressure)
     	.def("ComputeExactVelocity", &ManufacturedSolutionUtility::ComputeExactVelocity)
     	.def("ComputeExactPressure", &ManufacturedSolutionUtility::ComputeExactPressure)
     	.def("ComputeVelocityRelativeError", &ManufacturedSolutionUtility::ComputeVelocityRelativeError)
     	.def("ComputePressureRelativeError", &ManufacturedSolutionUtility::ComputePressureRelativeError)
     	;

}

} // namespace Python.
} // Namespace Kratos
