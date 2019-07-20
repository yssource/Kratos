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
#include "add_custom_manufactured_to_python.h"
#include "custom_manufactured/manufactured_solution.h"
#include "custom_manufactured/codina_vortex.h"
#include "custom_manufactured/eca_flow.h"


namespace Kratos {
namespace Python {

void  AddCustomManufacturedToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< ManufacturedSolution, ManufacturedSolution::Pointer>(m, "ManufacturedSolution")
        .def(py::init<Properties::Pointer&, Parameters::Pointer&>() )
        .def("BodyForce", &ManufacturedSolution::BodyForce)
        .def("Velocity", &ManufacturedSolution::Velocity)
        .def("Pressure", &ManufacturedSolution::Pressure)
        .def("TimeDerivative", &ManufacturedSolution::TimeDerivative)
        .def("ConvectiveTerm", &ManufacturedSolution::ConvectiveTerm)
        .def("ViscousTerm", &ManufacturedSolution::ViscousTerm)
        .def("PressureGradient", &ManufacturedSolution::PressureGradient)
        .def("VelocityGradient", &ManufacturedSolution::VelocityGradient)
        .def("VelocityLaplacian", &ManufacturedSolution::VelocityLaplacian)
        .def("Reynolds", &ManufacturedSolution::Reynolds)
        .def("Strouhal", &ManufacturedSolution::Strouhal)
        .def("GetParameters", &ManufacturedSolution::GetParameters)
        .def("GetProperties", &ManufacturedSolution::GetProperties)
        ;

    py::class_< CodinaVortex, CodinaVortex::Pointer, ManufacturedSolution>(m, "CodinaVortex")
        .def(py::init<Properties::Pointer&, Parameters::Pointer&>() )
        ;

    py::class_< EcaFlow, EcaFlow::Pointer, ManufacturedSolution>(m, "EcaFlow")
        .def(py::init<Properties::Pointer&, Parameters::Pointer&>() )
        ;

}

} // namespace Python.
} // Namespace Kratos
