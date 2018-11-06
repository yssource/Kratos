//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace pybind11;



  PYBIND11_MODULE(KratosCompressiblePotentialFlowApplication,m)
  {

	  class_<KratosCompressiblePotentialFlowApplication,
			  KratosCompressiblePotentialFlowApplication::Pointer,
			  KratosApplication >(m,"KratosCompressiblePotentialFlowApplication")
			  .def(init<>())
			;

    AddCustomProcessesToPython(m);
	//registering variables in python

//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,VELOCITY_INFINITY);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,VELOCITY_LOWER);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_LOWER);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,UPPER_SURFACE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,LOWER_SURFACE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,UPPER_WAKE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,LOWER_WAKE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,POTENTIAL_JUMP);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ENERGY_NORM_REFERENCE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,POTENTIAL_ENERGY_REFERENCE);

	


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
