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

#if !defined(KRATOS_ADD_CUSTOM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>
// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos
{
    namespace Python
    {
        void  AddCustomConstitutiveLawsToPython(pybind11::module& m);
    }  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED  defined
