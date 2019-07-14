//   __  __  _____ _____  _____                     _ _           _   _
//  |  \/  |/ ____|_   _|/ ____|  /\               | (_)         | | (_)
//  | \  / | |  __  | | | (___   /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |\/| | | |_ | | |  \___ \ / /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \
//  | |  | | |__| |_| |_ ____) / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//  |_|  |_|\_____|_____|_____/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                     | |   | |
//                                     |_|   |_|
//
//  License: BSD License
//   license: MGISApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

#if defined(KRATOS_PYTHON)
// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "mgis_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosMGISApplication, m)
{
    class_<KratosMGISApplication,
    KratosMGISApplication::Pointer,
    KratosApplication>(m, "KratosMGISApplication")
    .def(init<>());

    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomConstitutiveLawsToPython(m);

    // Registering variables in python
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, )
}


}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
