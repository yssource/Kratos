//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "manufactured_fluid_solutions_application.h"
#include "manufactured_fluid_solutions_application_variables.h"

namespace Kratos {

KratosManufacturedFluidSolutionsApplication::KratosManufacturedFluidSolutionsApplication():
    KratosApplication("ManufacturedFluidSolutionsApplication")
    {}

void KratosManufacturedFluidSolutionsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosManufacturedFluidSolutionsApplication..." << std::endl;

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EXACT_VELOCITY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_RELATIVE_ERROR )
    KRATOS_REGISTER_VARIABLE( EXACT_PRESSURE )
    KRATOS_REGISTER_VARIABLE( PRESSURE_RELATIVE_ERROR )
    KRATOS_REGISTER_VARIABLE( MANUFACTURED_PARAMETERS )

}

}  // namespace Kratos.
