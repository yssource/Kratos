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

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "rans_constitutive_laws_application.h"
#include "rans_constitutive_laws_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosRANSConstitutiveLawsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosRANSConstitutiveLawsApplication,
        KratosRANSConstitutiveLawsApplication::Pointer,
        KratosApplication>(m, "KratosRANSConstitutiveLawsApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SQRT_KINETIC_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SQRT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SQRT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_SQRT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_LOG_KINETIC_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_LOG_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_LOG_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_LOG_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_MODELLING_PROCESS_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_Y_PLUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTIAL_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_VELOCITY )


    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
