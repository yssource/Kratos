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

    typedef array_1d<double, 3> Vector3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssignMaterialOrientationUtility(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

    /// Destructor.
    virtual ~AssignMaterialOrientationUtility() = default;


    ///@}
    ///@name Operations
    ///@{

    void Execute(Parameters MethodParameters);

    void WriteFiberAngles(const std::string& rFileName = "MaterialOrientations.dat");

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

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    int mEchoLevel = 0;

    ///@}
    ///@name Private Operations
    ///@{

    void CheckAndReadVectors(Parameters ThisParameters, const std::string KeyName, Vector3& rVector);

    ///@}

}; // Class AssignMaterialOrientationUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ASSIGN_MATERIAL_ORIENTATION_UTILITY_H  defined
