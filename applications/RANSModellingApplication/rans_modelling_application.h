//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//


#if !defined(KRATOS_RANS_MODELLING_APPLICATION_H_INCLUDED )
#define  KRATOS_RANS_MODELLING_APPLICATION_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Primal elements
#include "custom_elements/evm_k_epsilon/evm_k_element.h"
#include "custom_elements/evm_k_epsilon/evm_epsilon_element.h"

// Adjoint elements
#include "custom_elements/evm_k_epsilon/evm_k_adjoint_element.h"
#include "custom_elements/evm_k_epsilon/evm_epsilon_adjoint_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_vms_adjoint_element.h"
#include "custom_elements/evm_k_epsilon/evm_monolithic_k_epsilon_vms_adjoint_element.h"

namespace Kratos {

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
class KratosRANSModellingApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosRANSModellingApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosRANSModellingApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosRANSModellingApplication();

    /// Destructor.
    ~KratosRANSModellingApplication() override {}

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
        return "KratosRANSModellingApplication";
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
          KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

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

    /// k-epsilon turbulence model adjoint elements
    const EvmKAdjointElement<2, 3> mRANSEVMKAdjoint2D;
    const EvmKAdjointElement<3, 4> mRANSEVMKAdjoint3D;

    const EvmEpsilonAdjointElement<2, 3> mRANSEVMEpsilonAdjoint2D;
    const EvmEpsilonAdjointElement<3, 4> mRANSEVMEpsilonAdjoint3D;

    const EvmKEpsilonVMSAdjointElement<2> mRANSEVMKEpsilonVMSAdjoint2D;
    const EvmKEpsilonVMSAdjointElement<3> mRANSEVMKEpsilonVMSAdjoint3D;

    const EvmMonolithicKEpsilonVMSAdjointElement<2> mRANSEVMMonolithicKEpsilonVMSAdjoint2D;
    const EvmMonolithicKEpsilonVMSAdjointElement<3> mRANSEVMMonolithicKEpsilonVMSAdjoint3D;

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
    KratosRANSModellingApplication& operator=(KratosRANSModellingApplication const& rOther);

    /// Copy constructor.
    KratosRANSModellingApplication(KratosRANSModellingApplication const& rOther);


    ///@}

}; // Class KratosRANSModellingApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_RANS_MODELLING_APPLICATION_H_INCLUDED  defined
