//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// #define KRATOS_CG_SOLVER_H_EXCLUDED

// System includes


// External includes
#include <pybind11/pybind11.h>



// Project includes
#include "includes/define_python.h"
#include "includes/kratos_version.h"
#include "add_testing_to_python.h"

namespace Kratos
{

namespace Python
{

char const* greet()
{
	std::stringstream header;
	header << "Hello, I am Kratos Multi-Physics Testing Module " << GetVersionString() <<" ;-)";
    return header.str().c_str();
}

PYBIND11_MODULE(KratosTests, m)
{
    namespace py = pybind11;

    AddTestingToPython(m);
}


}  // namespace Python.

}  // namespace Kratos.
