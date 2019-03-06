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
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"
#include "utilities/normal_calculation_utils.h"
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

    void InitializeConditionFlagsForModelPart(ModelPart* pModelPart, const Flags& rFlag)
    {
        KRATOS_TRY

        auto& conditions_array = pModelPart->Conditions();

        for (auto& it_cond : conditions_array)
        {
            bool is_flag_true = true;

            auto& r_geometry = it_cond.GetGeometry();

            for (auto& it_node : r_geometry.Points())
                if (!(it_node.Is(rFlag)))
                    is_flag_true = false;
            it_cond.Set(rFlag, is_flag_true);
        }

        KRATOS_CATCH("");
    }

    void FindConditionsParentElements(ModelPart* pModelPart)
    {
        FindNodalNeighboursProcess find_nodal_neighbours_process(*pModelPart);
        find_nodal_neighbours_process.Execute();

        const int number_of_conditions = pModelPart->NumberOfConditions();

#pragma omp parallel for
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
        {
            Condition& r_condition = *(pModelPart->ConditionsBegin() + i_cond);

            Geometry<Node<3>>& r_geometry = r_condition.GetGeometry();
            WeakPointerVector<Element> element_candidates;

            const IndexType number_of_nodes = r_geometry.PointsNumber();

            for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                WeakPointerVector<Element>& r_node_element_candidates =
                    r_geometry[i_node].GetValue(NEIGHBOUR_ELEMENTS);
                for (IndexType i_elem = 0; i_elem < r_node_element_candidates.size(); ++i_elem)
                    element_candidates.push_back(r_node_element_candidates(i_elem));
            }

            std::vector<IndexType> node_ids(number_of_nodes);
            for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
                node_ids[i_node] = r_geometry[i_node].Id();

            std::sort(node_ids.begin(), node_ids.end());

            std::vector<IndexType> element_node_ids;
            for (IndexType i_elem = 0; i_elem < element_candidates.size(); ++i_elem)
            {
                Geometry<Node<3>>& r_element_geometry =
                    element_candidates[i_elem].GetGeometry();
                const IndexType number_of_element_nodes =
                    r_element_geometry.PointsNumber();
                element_node_ids.resize(number_of_element_nodes);

                for (IndexType i_node = 0; i_node < number_of_element_nodes; ++i_node)
                    element_node_ids[i_node] = r_element_geometry[i_node].Id();

                std::sort(element_node_ids.begin(), element_node_ids.end());

                // If all the node in the condition is included in the element, then it is the parent element for that condition
                if (std::includes(element_node_ids.begin(), element_node_ids.end(),
                                  node_ids.begin(), node_ids.end()))
                {
                    r_condition.SetValue(PARENT_ELEMENT, element_candidates(i_elem));
                    break;
                }
            }
        }
    }

    void FindConditionGaussPointIndices(ModelPart* pModelPart)
    {
        const int number_of_conditions = pModelPart->NumberOfConditions();

#pragma omp parallel for
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
        {
            Condition& r_condition = *(pModelPart->ConditionsBegin() + i_cond);
            Geometry<Node<3>>& r_condition_geometry = r_condition.GetGeometry();
            const Element& r_element = *(r_condition.GetValue(PARENT_ELEMENT).lock());
            const Geometry<Node<3>>& r_element_geometry = r_element.GetGeometry();

            Vector gauss_weights;
            Matrix shape_functions;
            Geometry<Node<3>>::ShapeFunctionsGradientsType shape_derivatives;

            CalculateGeometryData(r_element_geometry, r_element.GetIntegrationMethod(),
                                  gauss_weights, shape_functions, shape_derivatives);

            const unsigned int num_gauss_points = gauss_weights.size();
            const unsigned int number_of_element_nodes =
                r_element_geometry.PointsNumber();
            const unsigned int number_of_condition_nodes =
                r_condition_geometry.PointsNumber();

            std::vector<int>& r_gauss_point_indices =
                r_condition.GetValue(GAUSS_POINT_INDICES);
            r_gauss_point_indices.clear();
            r_gauss_point_indices.resize(number_of_condition_nodes);

            for (unsigned int i_cond_node = 0;
                 i_cond_node < number_of_condition_nodes; ++i_cond_node)
            {
                int g_index = -1;
                double min_distance = 0.0;
                for (unsigned int g = 0; g < num_gauss_points; g++)
                {
                    const Matrix& r_shape_derivatives = shape_derivatives[g];
                    const Vector& gauss_shape_functions = row(shape_functions, g);

                    double gp_x(0.0), gp_y(0.0), gp_z(0.0);

                    for (unsigned int i_node = 0; i_node < number_of_element_nodes; ++i_node)
                    {
                        gp_x += gauss_shape_functions[i_node] *
                                r_element_geometry[i_node].X();
                        gp_y += gauss_shape_functions[i_node] *
                                r_element_geometry[i_node].Y();
                        gp_z += gauss_shape_functions[i_node] *
                                r_element_geometry[i_node].Z();
                    }

                    double distance = 0.0;
                    distance +=
                        std::pow(gp_x - r_condition_geometry[i_cond_node].X(), 2);
                    distance +=
                        std::pow(gp_y - r_condition_geometry[i_cond_node].Y(), 2);
                    distance +=
                        std::pow(gp_z - r_condition_geometry[i_cond_node].Z(), 2);
                    distance = std::sqrt(distance);

                    if (g_index == -1)
                    {
                        g_index = 0;
                        min_distance = distance;
                    }

                    if (min_distance > distance)
                    {
                        min_distance = distance;
                        g_index = g;
                    }
                }
                r_gauss_point_indices[i_cond_node] = g_index;
            }
        }
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

    virtual void InitializeConditions()
    {
        KRATOS_THROW_ERROR(std::runtime_error,
                           "Calling base class InitializeConditions.", "");
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

            y_plus = EvmKepsilonModelUtilities::CalculateYplus(
                velocity_norm, wall_distance, nu, von_karman, beta, 10);
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
