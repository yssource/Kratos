//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//                   Philipp Bucher
//


#if !defined(KRATOS_ASSIGN_MATERIAL_ORIENTATION_UTILITY_H )
#define  KRATOS_ASSIGN_MATERIAL_ORIENTATION_UTILITY_H


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
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

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AssignMaterialOrientationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignMaterialOrientationUtility
    KRATOS_CLASS_POINTER_DEFINITION(AssignMaterialOrientationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssignMaterialOrientationUtility(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

    /// Destructor.
    virtual ~AssignMaterialOrientationUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute(Parameters MethodParameters);

    void WriteFiberAngles(const std::string& rFileName = "MaterialOrientations.dat");

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
    virtual std::string Info() const {}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPart& mrModelPart;
    int mEchoLevel = 0;

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

    // /// Assignment operator.
    // AssignMaterialOrientationUtility& operator=(AssignMaterialOrientationUtility const& rOther);

    // /// Copy constructor.
    // AssignMaterialOrientationUtility(AssignMaterialOrientationUtility const& rOther);


    ///@}

}; // Class AssignMaterialOrientationUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                AssignMaterialOrientationUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const AssignMaterialOrientationUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ASSIGN_MATERIAL_ORIENTATION_UTILITY_H  defined
