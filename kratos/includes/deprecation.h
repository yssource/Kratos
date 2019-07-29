#ifndef KRATOS_DEPRECATION_H_INCLUDED
#define KRATOS_DEPRECATION_H_INCLUDED

#include "define.h"
#include "includes/kratos_version.h"

namespace Kratos {

template<unsigned int Major, unsigned int Minor, unsigned int Patch>
struct CurrentVersion
{
    constexpr static unsigned int major = Major;
    constexpr static unsigned int minor = Minor;
    constexpr static unsigned int patch = Patch;
};

template<unsigned int Major, unsigned int Minor, unsigned int Patch>
struct DeprecatedInVersion
{
    constexpr static unsigned int major = Major;
    constexpr static unsigned int minor = Minor;
    constexpr static unsigned int patch = Patch;
};

template<class A, class B>
struct IsNewerThan
{
    constexpr static bool value = (A::major > B::major) || ((A::major == B::major) && ((A::minor > B::minor) || ( (A::minor == B::minor) && (A::patch >= B::patch) )));
};

template<class CurrentVersion, class DeprecatedSinceVersion, bool Deprecated = IsNewerThan<CurrentVersion,DeprecatedSinceVersion>::value >
struct CheckDeprecatedCall {};

template<class CurrentVersion, class DeprecatedSinceVersion>
struct CheckDeprecatedCall<CurrentVersion,DeprecatedSinceVersion,false>
{
    void AccessingDeprecatedFunction();

    template<class DeprecatedSince>
    KRATOS_DEPRECATED static void AccessingDeprecatedFunction()
    {
        //#pragma message("Warning: Calling function that will be deprecated.")
    }
};

template<class CurrentVersion, class DeprecationVersion>
void AccessingDeprecatedFunction()
{
    CheckDeprecatedCall<CurrentVersion, DeprecationVersion>::template AccessingDeprecatedFunction<DeprecationVersion>();
}

#ifndef KRATOS_DEPRECATED_FROM_VERSION_INTERNAL
#define KRATOS_DEPRECATED_FROM_VERSION_INTERNAL(current_major, current_minor, current_patch, version_major, version_minor, version_patch) \
{ \
using CurrentKratosVersion = CurrentVersion<current_major, current_minor, current_patch>;          \
using DeprecationVersion = DeprecatedInVersion<version_major, version_minor, version_patch>;          \
CheckDeprecatedCall<CurrentKratosVersion, DeprecationVersion>::template AccessingDeprecatedFunction<DeprecationVersion>(); \
}    \

#endif


#ifndef KRATOS_DEPRECATED_FROM_VERSION
#define KRATOS_DEPRECATED_FROM_VERSION(version_major, version_minor, version_patch) \
KRATOS_DEPRECATED_FROM_VERSION_INTERNAL(KRATOS_MAJOR_VERSION, KRATOS_MINOR_VERSION, 0, version_major, version_minor, version_patch)
#endif

}

#endif