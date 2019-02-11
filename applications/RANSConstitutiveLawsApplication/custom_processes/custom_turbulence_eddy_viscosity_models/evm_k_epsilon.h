//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					     Kratos default license: kratos/license.txt
//
//  Main author: Suneth Warnakulasuriya
//

#if !defined(KRATOS_TURBULENCE_EVM_K_EPSILON_H_INCLUDED)
#define KRATOS_TURBULENCE_EVM_K_EPSILON_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "eddy_viscosity_model.h"



namespace Kratos
{
///@addtogroup RANSConstitutiveLawsApplication
///@{

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

/// Auxiliary process to set Boussinesq buoyancy forces in variable temperature flows.
/** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis
    so that the fluid element can take natural convection into account.

    This process makes use of the following data:
    - TEMPERATURE from the nodal solution step data: current temperature for the node (mandatory).
    - AMBIENT_TEMPERATURE from ProcessInfo: The reference temperature for the simulation (mandatory).
    - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
    - thermal_expansion_coefficient from the Parameters: a double defining the thermal expansion coefficient for the fluid (optional).

    With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
    The force is set to (1 + thermal_expansion_coefficient*(temperature - ambient_temperature) ) * g

    If the thermal expansion coefficient is not provided, it is assumed to be (1/ambient_temperature).
    This is the usual value for perfect gases (if the temperature is given in Kelvin).
 */

class TurbulenceEvmKepsilon : public TurbulenceEddyViscosityModel
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TurbulenceEvmKepsilon
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEvmKepsilon);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TurbulenceEvmKepsilon(ModelPart& rModelPart,
                          Parameters& rParameters,
                          std::vector<ModelPart>& rWallConditionsModelPartList,
                          std::vector<ModelPart>& rInletConditionsModelPartList,
                          std::vector<ModelPart>& rOutletConditionsModelPartList)
        : TurbulenceEddyViscosityModel(rModelPart, rParameters, rWallConditionsModelPartList, rInletConditionsModelPartList, rOutletConditionsModelPartList)
    {
    }

    /// Destructor.
    ~TurbulenceEvmKepsilon() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual Element const& GetReferenceElement()
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "Calling turbulence modelling base class "
                           "GetReferenceElement method",
                           "");

        return Element();

        KRATOS_CATCH("");
    }

    virtual Condition const& GetReferenceCondition()
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "Calling turbulence modelling base class "
                           "GetReferenceCondition method",
                           "");

        return Condition();

        KRATOS_CATCH("");
    }

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    ///@}
    ///@name Member Variables
    ///@{

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
    TurbulenceEvmKepsilon& operator=(TurbulenceEvmKepsilon const& rOther);

    /// Copy constructor.
    TurbulenceEvmKepsilon(TurbulenceEvmKepsilon const& rOther);

    ///@}

}; // Class TurbulenceEvmKepsilon

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(std::ostream& rOStream, const TurbulenceEvmKepsilon& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TURBULENCE_EVM_K_EPSILON_H_INCLUDED  defined
