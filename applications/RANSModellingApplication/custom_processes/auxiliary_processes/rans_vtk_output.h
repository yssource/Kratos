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

#if !defined(KRATOS_RANS_VTK_OUTPUT_PROCESS_H_INCLUDED)
#define KRATOS_RANS_VTK_OUTPUT_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "input_output/vtk_output.h"

#include "custom_utilities/rans_check_utilities.h"
#include "custom_utilities/rans_variable_utils.h"

namespace Kratos
{
///@addtogroup RANSModellingApplication
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

/**
 * @brief Checks lower and upper bounds of a scalar
 *
 * This process checks lower and upper bounds of a variable in nodes in a given modelpart
 *
 */

class RansVtkOutputProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansVtkOutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansVtkOutputProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansVtkOutputProcess(Model& rModel, Parameters& rParameters)
        : Process(), mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "folder_name"                                 : "auxiliary_output"
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();
        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
        mpVtkOutput =  Kratos::shared_ptr<VtkOutput>(new VtkOutput(r_model_part, mrParameters));

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansVtkOutputProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        KRATOS_TRY

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mModelPartName);

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        mCount++;
        const int current_step = r_model_part.GetProcessInfo()[STEP];
        r_model_part.GetProcessInfo()[STEP] = mCount;

        mpVtkOutput->PrintOutput();

        r_model_part.GetProcessInfo()[STEP] = current_step;

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
    std::string Info() const override
    {
        return std::string("RansVtkOutputProcess");
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
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

    Model& mrModel;
    Parameters& mrParameters;
    std::string mModelPartName;
    long unsigned int mCount = 0;

    VtkOutput::Pointer mpVtkOutput;

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
    RansVtkOutputProcess& operator=(RansVtkOutputProcess const& rOther);

    /// Copy constructor.
    RansVtkOutputProcess(RansVtkOutputProcess const& rOther);

    ///@}

}; // Class RansVtkOutputProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansVtkOutputProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_VTK_OUTPUT_PROCESS_H_INCLUDED defined
