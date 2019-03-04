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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"
#include "rans_constitutive_laws_application_variables.h"

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

    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Pointer definition of TurbulenceEddyViscosityModelProcess
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEddyViscosityModelProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TurbulenceEddyViscosityModelProcess(ModelPart& rModelPart,
                                        Parameters& rParameters,
                                        typename TLinearSolver::Pointer pLinearSolver);

    /// Destructor.
    ~TurbulenceEddyViscosityModelProcess() override
    {
        // delete mpDistanceCalculator;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void Execute() override;

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
    ModelPart& mrModelPart;
    Parameters& mrParameters;
    typename TLinearSolver::Pointer mpLinearSolver;

    bool mIsMeshMoving;

    unsigned int mEchoLevel;
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void AddSolutionStepVariables()
    {
        KRATOS_TRY

        mrModelPart.GetNodalSolutionStepVariablesList().push_back(DISTANCE);
        mrModelPart.GetNodalSolutionStepVariablesList().push_back(FLAG_VARIABLE);
        mrModelPart.GetNodalSolutionStepVariablesList().push_back(KINEMATIC_VISCOSITY);
        mrModelPart.GetNodalSolutionStepVariablesList().push_back(TURBULENT_VISCOSITY);

        KRATOS_INFO("TurbulenceModel") << "Added eddy viscosity turbulence "
                                          "model solution step variables.\n";

        KRATOS_CATCH("");
    }

    virtual void AddDofs()
    {
        KRATOS_INFO("TurbulenceModel")
            << "Added eddy viscosity turbulence model dofs.\n";
    }

    virtual void InitializeTurbulenceModelPart()
    {
        KRATOS_THROW_ERROR(std::runtime_error,
                           "Calling base class InitializeTurbulenceModelPart.",
                           "");
    }

    virtual void UpdateFluidViscosity()
    {
        KRATOS_TRY

        NodesArrayType& nodes = mrModelPart.Nodes();

// Modifying viscosity of the nodes with the calculated turbulent viscosity
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            auto it_node = nodes.begin() + i;
            const double kinematic_viscosity =
                it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double turbulent_viscosity =
                it_node->FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            double& effective_viscosity = it_node->FastGetSolutionStepValue(VISCOSITY);
            effective_viscosity = kinematic_viscosity + turbulent_viscosity;
        }

        KRATOS_CATCH("");
    }

    virtual void InitializeConditionFlags(const Flags& rFlag)
    {
        KRATOS_THROW_ERROR(std::runtime_error,
                           "Calling base class InitializeConditionFlags.", "");
    }

    void GenerateModelPart(ModelPart& rOriginModelPart,
                           ModelPart& rDestinationModelPart,
                           const Element& rReferenceElement,
                           const Condition& rReferenceCondition);

    void CalculateYplus(unsigned int Step = 0)
    {
        int number_of_nodes = mrModelPart.NumberOfNodes();
        const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

        const double von_karman = r_current_process_info[WALL_VON_KARMAN];
        const double beta = r_current_process_info[WALL_SMOOTHNESS_BETA];

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; ++i)
        {
            Node<3>& r_node = *(mrModelPart.NodesBegin() + i);

            const array_1d<double, 3>& r_velocity =
                r_node.FastGetSolutionStepValue(VELOCITY, Step);
            const double velocity_norm = norm_2(r_velocity);

            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);

            double& y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);

            // try linear law
            y_plus = std::sqrt(velocity_norm * wall_distance / nu);

            // If the linear low doesnt match within the range, try logrithmic law
            if (y_plus > 11.06)
            {
                unsigned int max_u_tau_iterations = 10, i;
                double u_tau = std::sqrt(velocity_norm * nu / wall_distance);
                double prev_u_tau = 0.0;
                for (i = 0; i < max_u_tau_iterations; ++i)
                {
                    prev_u_tau = u_tau;
                    u_tau = velocity_norm /
                            (std::log(u_tau * wall_distance / nu) / von_karman + beta);
                }
                const double delta_u_tau = std::abs(u_tau - prev_u_tau);
                KRATOS_INFO_IF("TurbulenceEvmProcess", delta_u_tau > 1e-5) << "WARNING: Maximum number of iterations reached for y_plus calculation. error_u_tau = "
                                                                           << std::scientific
                                                                           << delta_u_tau
                                                                           << ".\n";
                y_plus = u_tau * wall_distance / nu;
            }
        }
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

    VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>* mpDistanceCalculator;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances();

    void InitializeNodeFlags(const Parameters& rParameters, const Flags& rFlag)
    {
        KRATOS_TRY

        for (std::size_t i = 0; i < rParameters.size(); ++i)
        {
            std::string model_part_name = rParameters.GetArrayItem(i).GetString();
            KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(model_part_name))
                << "TurbulenceEddyViscosityModelProcess: Wall condition "
                << model_part_name << " not found." << std::endl;
            ModelPart& current_model_part = mrModelPart.GetSubModelPart(model_part_name);

            NodesArrayType& nodes_array = current_model_part.Nodes();

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
            {
                auto it_node = nodes_array.begin() + i;
                it_node->Set(rFlag, true);
            }
        }

        KRATOS_CATCH("");
    }

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
