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
        .def("ComputeExactMaterialAcceleration", &ManufacturedSolutionUtility::ComputeExactMaterialAcceleration)
     	.def("ComputeVelocityRelativeError", &ManufacturedSolutionUtility::ComputeVelocityRelativeError)
     	.def("ComputePressureRelativeError", &ManufacturedSolutionUtility::ComputePressureRelativeError)
        .def("ComputeMaterialAccelerationError", &ManufacturedSolutionUtility::ComputeMaterialAccelerationError)
        .def("ComputeError", &ManufacturedSolutionUtility::ComputeError<Variable<double>>)
        .def("ComputeError", &ManufacturedSolutionUtility::ComputeError<Variable<array_1d<double,3>>>)
        .def("ComputeError", &ManufacturedSolutionUtility::ComputeError<Variable<array_1d<double,4>>>)
        .def("ComputeError", &ManufacturedSolutionUtility::ComputeError<Variable<array_1d<double,6>>>)
        .def("ComputeError", &ManufacturedSolutionUtility::ComputeError<Variable<array_1d<double,9>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<Variable<double>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<Variable<array_1d<double,3>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<Variable<array_1d<double,4>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<Variable<array_1d<double,6>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<Variable<array_1d<double,9>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("ComputeMean", &ManufacturedSolutionUtility::ComputeMean<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<Variable<double>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<Variable<array_1d<double,3>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<Variable<array_1d<double,4>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<Variable<array_1d<double,6>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<Variable<array_1d<double,9>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("ComputeRootMeanSquare", &ManufacturedSolutionUtility::ComputeRootMeanSquare<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("RecoverMaterialAcceleration", &ManufacturedSolutionUtility::RecoverMaterialAcceleration)
        .def("BDF1", &ManufacturedSolutionUtility::BDF1<Variable<double>>)
        .def("BDF1", &ManufacturedSolutionUtility::BDF1<Variable<array_1d<double,3>>>)
        .def("BDF1", &ManufacturedSolutionUtility::BDF1<Variable<array_1d<double,4>>>)
        .def("BDF1", &ManufacturedSolutionUtility::BDF1<Variable<array_1d<double,6>>>)
        .def("BDF1", &ManufacturedSolutionUtility::BDF1<Variable<array_1d<double,9>>>)
        .def("ComputeNodalCFL", &ManufacturedSolutionUtility::ComputeNodalCFL)
     	;

}

} // namespace Python.
} // Namespace Kratos
