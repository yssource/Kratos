//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#include "manufactured_fluid_solutions_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EXACT_VELOCITY )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_RELATIVE_ERROR )
KRATOS_CREATE_VARIABLE( double, EXACT_PRESSURE )
KRATOS_CREATE_VARIABLE( double, PRESSURE_RELATIVE_ERROR )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EXACT_CONVECTION )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONVECTION_ERROR )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EXACT_MATERIAL_ACCELERATION )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ACCELERATION_ERROR )
}
