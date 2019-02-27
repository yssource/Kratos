//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Suneth Warnakulasuriya
//
//

#if !defined(KRATOS_ELEMENT_DERIVATIVES_EXTENSION_INCLUDED)
#define KRATOS_ELEMENT_DERIVATIVES_EXTENSION_INCLUDED

// System includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class ElementDerivativesExtension
 * @ingroup KratosCore
 * @brief Interface extensions for elements and conditions.
 */
class ElementDerivativesExtension
{
public:
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    KRATOS_CLASS_POINTER_DEFINITION(ElementDerivativesExtension);

    virtual ~ElementDerivativesExtension()
    {
    }

    virtual void GetFirstDerivativesDofList(DofsVectorType& rElementalDofList,
                                            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "ElementDerivativesExtension: Calling base class GetFirstDerivativesDofList method.", "");
    }

    virtual void GetSecondDerivativesDofList(DofsVectorType& rElementalDofList,
                                             ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "ElementDerivativesExtension: Calling base class GetSecondDerivativesDofList method.", "");
    }

    virtual std::ostream& Print(std::ostream& os) const
    {
        return os;
    }

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }
};

///@} // Kratos Classes

} // namespace Kratos.

#endif // KRATOS_ELEMENT_DERIVATIVES_EXTENSION_INCLUDED  defined
