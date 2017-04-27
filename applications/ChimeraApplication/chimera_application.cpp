//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "chimera_application.h"
#include "chimera_application_variables.h"


namespace Kratos {

KratosChimeraApplication::KratosChimeraApplication() {}

void KratosChimeraApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosChimeraApplication... " << std::endl;

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CHIM_NEUMANN_COND )
   KRATOS_REGISTER_VARIABLE(MPC_POINTER);


}
}  // namespace Kratos.
