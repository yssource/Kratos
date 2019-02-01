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

}
}  // namespace Kratos.
