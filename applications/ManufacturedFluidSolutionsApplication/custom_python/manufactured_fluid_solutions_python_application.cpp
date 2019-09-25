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

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "manufactured_fluid_solutions_application.h"
#include "manufactured_fluid_solutions_application_variables.h"
#include "custom_python/add_custom_manufactured_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosManufacturedFluidSolutionsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosManufacturedFluidSolutionsApplication,
        KratosManufacturedFluidSolutionsApplication::Pointer,
        KratosApplication>(m, "KratosManufacturedFluidSolutionsApplication")
        .def(py::init<>())
        ;

    AddCustomManufacturedToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXACT_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_RELATIVE_ERROR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXACT_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_RELATIVE_ERROR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXACT_CONVECTION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONVECTION_ERROR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXACT_MATERIAL_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MATERIAL_ACCELERATION_ERROR )
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
