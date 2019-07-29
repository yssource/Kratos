
#include "testing/testing.h"
#include "includes/deprecation.h"
#include "includes/kratos_version.h"

namespace Kratos {
namespace Testing {

void DeprecatedFunction()
{
    // deprecated from current version
    KRATOS_DEPRECATED_FROM_VERSION(GetMajorVersion(),GetMinorVersion(),0);
    // deprecated in next major release
    //KRATOS_DEPRECATED_FROM_VERSION(GetMajorVersion()+1,GetMinorVersion(),0);
}

}
}