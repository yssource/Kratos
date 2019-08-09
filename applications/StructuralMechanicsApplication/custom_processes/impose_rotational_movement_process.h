// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_IMPOSE_ROTATIONAL_MOVEMENT_PROCESS)
#define KRATOS_IMPOSE_ROTATIONAL_MOVEMENT_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/python_function_callback_utility.h"

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

/**
 * @class ImposeRotationalMovementProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This method assign linear kinematic constrains to a certain submodelpart to impose a rotational movement
 * @details It assigns to the given point or it detects the center of gravity automatically
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ImposeRotationalMovementProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImposeRotationalMovementProcess
    KRATOS_CLASS_POINTER_DEFINITION(ImposeRotationalMovementProcess);

    /// General type definitions
    typedef Node<3>                                                      NodeType;

    /// General containers type definitions
    typedef ModelPart::MasterSlaveConstraintContainerType ConstraintContainerType;

    /// The DoF type definition
    typedef Dof<double>                                                   DofType;

    /// The DoF pointer vector type definition
    typedef std::vector< DofType::Pointer >                  DofPointerVectorType;

    /// Auxiliar matrix and vector
    typedef Matrix                                                     MatrixType;
    typedef Vector                                                     VectorType;

    /// Definitions of the integers
    typedef std::size_t                                                 IndexType;
    typedef std::size_t                                                  SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rThisModelPart The model part to compute
     * @param ThisParameters The parameters of configuration
     */
    ImposeRotationalMovementProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ImposeRotationalMovementProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

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
        return "ImposeRotationalMovementProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ImposeRotationalMovementProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrThisModelPart;                       /// The model part to compute
    Parameters mThisParameters;                       /// The parameters (can be used for general pourposes)
    NodeType::Pointer mpMasterNode;                   /// The master node which is considered to impose the rotational movement
    PythonGenericFunctionUtility::Pointer mpFunction; /// The python function used, depends on X, Y, Z, and t

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
    ImposeRotationalMovementProcess& operator=(ImposeRotationalMovementProcess const& rOther) = delete;

    /// Copy constructor.
    //ImposeRotationalMovementProcess(ImposeRotationalMovementProcess const& rOther);


    ///@}

}; // Class ImposeRotationalMovementProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ImposeRotationalMovementProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ImposeRotationalMovementProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}
#endif /* KRATOS_IMPOSE_ROTATIONAL_MOVEMENT_PROCESS defined */
