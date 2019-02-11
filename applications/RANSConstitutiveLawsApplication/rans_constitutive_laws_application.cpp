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
#include "rans_constitutive_laws_application.h"
#include "rans_constitutive_laws_application_variables.h"


namespace Kratos {

KratosRANSConstitutiveLawsApplication::KratosRANSConstitutiveLawsApplication():
    KratosApplication("RANSConstitutiveLawsApplication")
    {}

void KratosRANSConstitutiveLawsApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosRANSConstitutiveLawsApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_2 )

    // Turbulence model constants
    KRATOS_REGISTER_VARIABLE( WALL_SMOOTHNESS_BETA )
    KRATOS_REGISTER_VARIABLE( WALL_VON_KARMAN )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C_MU )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C1 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C2 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_MIXING_LENGTH )
    KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY_FRACTION )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )

    KRATOS_REGISTER_VARIABLE( ELEMENT_DERIVATIVES_DOFS_EXTENSION )
}
}  // namespace Kratos.
