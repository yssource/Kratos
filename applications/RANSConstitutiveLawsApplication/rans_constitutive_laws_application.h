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

#if !defined(KRATOS_RANS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED)
#define KRATOS_RANS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/evm_k_epsilon_wall_condition.h"
#include "custom_elements/evm_k_epsilon/evm_k_element.h"
#include "custom_elements/evm_k_epsilon/evm_epsilon_element.h"
#include "includes/kratos_application.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KratosRANSConstitutiveLawsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosRANSConstitutiveLawsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosRANSConstitutiveLawsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosRANSConstitutiveLawsApplication();

    /// Destructor.
    ~KratosRANSConstitutiveLawsApplication() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosRANSConstitutiveLawsApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    /// k-epsilon turbulence model elements
    const EvmKElement<2, 3> mRANSEVMK2D;
    const EvmKElement<3, 4> mRANSEVMK3D;

    const EvmEpsilonElement<2, 3> mRANSEVMEPSILON2D;
    const EvmEpsilonElement<3, 4> mRANSEVMEPSILON3D;

    /// k-epsilon turbulence model conditions
    const EvmKEpsilonWallCondition<2, 2> mKEpsilonWallCondition2D;
    const EvmKEpsilonWallCondition<3, 3> mKEpsilonWallCondition3D;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosRANSConstitutiveLawsApplication& operator=(KratosRANSConstitutiveLawsApplication const& rOther);

    /// Copy constructor.
    KratosRANSConstitutiveLawsApplication(KratosRANSConstitutiveLawsApplication const& rOther);

    ///@}

}; // Class KratosRANSConstitutiveLawsApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED  defined
