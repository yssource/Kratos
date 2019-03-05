//
//   Project Name:        KratosPlasmaDynamicsApplication $
//   Last Modified by:    $Author:    Marc Chung To Sang $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "plasma_dynamics_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosPlasmaDynamicsApplication, m)
{
    py::class_<KratosPlasmaDynamicsApplication,
    KratosPlasmaDynamicsApplication::Pointer,
    KratosApplication>(m, "KratosPlasmaDynamicsApplication")
    .def( py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomUtilitiesToPython(m);

    //Registering variables in python
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
