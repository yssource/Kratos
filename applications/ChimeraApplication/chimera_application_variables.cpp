#include "chimera_application_variables.h"

namespace Kratos
{
typedef MpcData::Pointer MpcDataPointerType;

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CHIM_NEUMANN_COND );

KRATOS_CREATE_VARIABLE(MpcDataPointerType, MPC_POINTER) // Amap of the master nodes to their corresponding weights
}
