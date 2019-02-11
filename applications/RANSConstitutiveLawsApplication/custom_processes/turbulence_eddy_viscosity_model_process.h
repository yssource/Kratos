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

#if !defined(KRATOS_TURBULENCE_EDDY_VISCOSITY_PROCESS_H_INCLUDED)
#define KRATOS_TURBULENCE_EDDY_VISCOSITY_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "../rans_constitutive_laws_application_variables.h"
#include "custom_turbulence_eddy_viscosity_models/eddy_viscosity_model.h"
#include "includes/define.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"
#include "includes/cfd_variables.h"
#include "../custom_strategies/general_convergence_criteria.h"
#include "../custom_strategies/residual_based_bossak_velocity_scheme.h"

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
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class TurbulenceEddyViscosityModelProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> ModelType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Pointer definition of TurbulenceEddyViscosityModelProcess
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEddyViscosityModelProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TurbulenceEddyViscosityModelProcess(ModelPart& rModelPart,
                                        Parameters& rParameters,
                                        TLinearSolver& rLinearSolver);

    /// Destructor.
    ~TurbulenceEddyViscosityModelProcess() override
    {
        delete mpEddyViscosityModel;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override
    {
      mpEddyViscosityModel->ExecuteBeforeSolutionLoop();
    }

    void ExecuteInitializeSolutionStep() override;

    void Execute() override;

    void ExecuteFinalizeSolutionStep() override
    {
      mpEddyViscosityModel->ExecuteFinalizeSolutionStep();
    }

    void ExecuteBeforeOutputStep() override
    {
      mpEddyViscosityModel->ExecuteBeforeOutputStep();
    }

    void ExecuteAfterOutputStep() override
    {
      mpEddyViscosityModel->ExecuteAfterOutputStep();
    }

    void ExecuteFinalize()
    {
      mpEddyViscosityModel->ExecuteFinalize();
    }

    int Check() override
    {
      return mpEddyViscosityModel->Check();
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

    void InitializeTurbulenceModelPart()
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(mrModelPart.HasSubModelPart("TurbulenceModelPart")) << "TurbulenceEddyViscosityModelProcess: TurbulenceModelPart is already found."
                                                                            << std::endl;
        mrModelPart.CreateSubModelPart("TurbulenceModelPart");
        mrTurbulenceModelPart = mrModelPart.GetSubModelPart("TurbulenceModelPart");

        const Element& rElem = mpEddyViscosityModel->GetReferenceElement();
        const Condition& rCond = mpEddyViscosityModel->GetReferenceCondition();

        this->GenerateModelPart(mrModelPart, mrTurbulenceModelPart, rElem, rCond);

        KRATOS_CATCH("");
    }

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
    Parameters& mrParameters;
    TLinearSolver& mrLinearSolver;
    TurbulenceEddyViscosityModel* mpEddyViscosityModel;

    ModelPart mrTurbulenceModelPart;

    VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> mDistanceCalculator;

    std::vector<ModelPart> mInletConditionsList;
    std::vector<ModelPart> mOutletConditionsList;
    std::vector<ModelPart> mWallConditionsList;

    bool mIsMeshMoving;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances();

    // void AssignBoundaryConditions();

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
    TurbulenceEddyViscosityModelProcess& operator=(
        TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    TurbulenceEddyViscosityModelProcess(
        TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    ///@}

}; // Class TurbulenceEddyViscosityModelProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TURBULENCE_EDDY_VISCOSITY_PROCESS_H_INCLUDED  defined
