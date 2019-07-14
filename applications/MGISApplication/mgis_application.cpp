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


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "mgis_application.h"

namespace Kratos
{

KratosMGISApplication::KratosMGISApplication():
    KratosApplication("MGISApplication")
    {}

void KratosMGISApplication::Register()
{
    // Calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosMGISApplication... " << std::endl;
}

}  // namespace Kratos.
