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

#if !defined(KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"
#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{
    typedef MpcData::Pointer MpcDataPointerType;
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( CHIM_NEUMANN_COND );
KRATOS_DEFINE_VARIABLE(MpcDataPointerType, MPC_POINTER);


}

#endif	/* KRATOS_CHIMERA_APPLICATION_VARIABLES_H_INCLUDED */
