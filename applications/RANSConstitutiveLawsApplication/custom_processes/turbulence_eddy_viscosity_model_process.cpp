#include "turbulence_eddy_viscosity_model_process.h"

#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
/* Public functions *******************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::TurbulenceEddyViscosityModelProcess(
    ModelPart& rModelPart, Parameters& rParameters, typename TLinearSolver::Pointer pLinearSolver)
    : mrModelPart(rModelPart), mrParameters(rParameters), mpLinearSolver(pLinearSolver)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "inlet_conditions"      : ["PLEASE_SPECIFY_INLET_CONDITIONS"],
        "outlet_conditions"     : ["PLEASE_SPECIFY_OUTLET_CONDITIONS"],
        "wall_conditions"       : ["PLEASE_SPECIFY_WALL_CONDITIONS"],
        "max_distance_calculation_iterations" : 2,
        "mesh_moving"       : false,
        "echo_level"        : 0,
        "model_properties"  : {}
    })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mIsMeshMoving = mrParameters["mesh_moving"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();

    KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitialize()
{
    this->InitializeTurbulenceModelPart();

    this->InitializeNodeFlags(mrParameters["inlet_conditions"], INLET);
    this->InitializeNodeFlags(mrParameters["outlet_conditions"], OUTLET);
    this->InitializeNodeFlags(mrParameters["wall_conditions"], STRUCTURE);
    this->InitializeNodeFlags(mrParameters["wall_conditions"], INLET, false);
    this->InitializeNodeFlags(mrParameters["wall_conditions"], OUTLET, false);

    NormalCalculationUtils normal_calculation_utility;
    normal_calculation_utility.CalculateOnSimplex(mrModelPart, TDim);

    this->InitializeConditions();

    mpDistanceCalculator =
        new VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
            mrModelPart, mpLinearSolver,
            mrParameters["max_distance_calculation_iterations"].GetInt());

    // Calculate the distances only once if the mesh is not moving
    CalculateWallDistances();

    NodesArrayType& nodes = mrModelPart.Nodes();

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
    {
        auto it_node = nodes.begin() + i;
        const double kinematic_viscosity = it_node->FastGetSolutionStepValue(VISCOSITY);
        it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = kinematic_viscosity;
    }
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    // Calculate the wall distances if the mesh has moved
    if (mIsMeshMoving)
        CalculateWallDistances();

    KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Execute()
{
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
std::string TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Info() const
{
    return "TurbulenceEddyViscosityModelProcess";
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "TurbulenceEddyViscosityModelProcess";
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PrintData(
    std::ostream& rOStream) const
{
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::CalculateWallDistances()
{
    KRATOS_TRY

    if (mEchoLevel > 0)
        KRATOS_INFO("TurbulenceModel") << "Calculating wall distances..." << std::endl;

    // Fixing the wall boundaries for wall distance calculation
    NodesArrayType& nodes_array = this->mrModelPart.Nodes();

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
    {
        auto it_node = nodes_array.begin() + i;
        // KRATOS_INFO("TurbulenceModel")<<it_node->Id()<<", "<<it_node->Is(STRUCTURE) <<std::endl;
        if (it_node->Is(STRUCTURE))
            it_node->FastGetSolutionStepValue(DISTANCE) = 0.0;
        else
            it_node->FastGetSolutionStepValue(DISTANCE) = 1.0;
    }

    mpDistanceCalculator->Execute();

    if (mEchoLevel > 0)
        KRATOS_INFO("TurbulenceModel")
            << "Finished calculating wall distances." << std::endl;

    KRATOS_CATCH("");
}

/* Protected functions ****************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::GenerateModelPart(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Element& rReferenceElement,
    const Condition& rReferenceCondition)
{
    // KRATOS_TRY

    // Copy general ModelPart properites
    // rDestinationModelPart.GetNodalSolutionStepVariablesList() =
    //     rOriginModelPart.GetNodalSolutionStepVariablesList();

    // rDestinationModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());
    rDestinationModelPart.SetProcessInfo(rOriginModelPart.pGetProcessInfo());
    rDestinationModelPart.SetProperties(rOriginModelPart.pProperties());

    // Copy tables
    rDestinationModelPart.Tables() = rOriginModelPart.Tables();

    // Copy the node list so that both model parts share the same nodes
    rDestinationModelPart.SetNodes(rOriginModelPart.pNodes());

    /* Create a new communicator for rDestinationModelPart and fill it with the information of the original one
     * Only "general" information and node lists are passed, element and condition lists will be created later
     * using the new elements.
     */

    Communicator& rReferenceComm = rOriginModelPart.GetCommunicator();
    Communicator::Pointer pDestinationComm = rReferenceComm.Create();
    pDestinationComm->SetNumberOfColors(rReferenceComm.GetNumberOfColors());
    pDestinationComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();
    pDestinationComm->LocalMesh().SetNodes(rReferenceComm.LocalMesh().pNodes());
    pDestinationComm->InterfaceMesh().SetNodes(rReferenceComm.InterfaceMesh().pNodes());
    pDestinationComm->GhostMesh().SetNodes(rReferenceComm.GhostMesh().pNodes());
    for (unsigned int i = 0; i < rReferenceComm.GetNumberOfColors(); i++)
    {
        pDestinationComm->pLocalMesh(i)->SetNodes(rReferenceComm.pLocalMesh(i)->pNodes());
        pDestinationComm->pInterfaceMesh(i)->SetNodes(
            rReferenceComm.pInterfaceMesh(i)->pNodes());
        pDestinationComm->pGhostMesh(i)->SetNodes(rReferenceComm.pGhostMesh(i)->pNodes());
    }

    rDestinationModelPart.SetCommunicator(pDestinationComm);

    // Reset element container and create new elements
    rDestinationModelPart.Elements().clear();
    rDestinationModelPart.Elements().reserve(rOriginModelPart.NumberOfElements());

    for (ModelPart::ElementsContainerType::iterator iEl = rOriginModelPart.ElementsBegin();
         iEl != rOriginModelPart.ElementsEnd(); iEl++)
    {
        Properties::Pointer pProp = iEl->pGetProperties();
        Element::Pointer pElem =
            rReferenceElement.Create(iEl->Id(), iEl->GetGeometry(), pProp);
        rDestinationModelPart.Elements().push_back(pElem);
    }

    // All elements are passed as local elements to the new communicator
    ModelPart::ElementsContainerType& rDestinationLocalElements =
        pDestinationComm->LocalMesh().Elements();
    rDestinationLocalElements.clear();
    rDestinationLocalElements.reserve(rDestinationModelPart.NumberOfElements());
    for (ModelPart::ElementsContainerType::ptr_iterator iEl =
             rDestinationModelPart.Elements().ptr_begin();
         iEl != rDestinationModelPart.Elements().ptr_end(); iEl++)
    {
        rDestinationLocalElements.push_back(*iEl);
    }

    // Reset condition container and create new conditions
    rDestinationModelPart.Conditions().clear();
    rDestinationModelPart.Conditions().reserve(rOriginModelPart.NumberOfConditions());
    for (ModelPart::ConditionsContainerType::iterator iCo = rOriginModelPart.ConditionsBegin();
         iCo != rOriginModelPart.ConditionsEnd(); iCo++)
    {
        Properties::Pointer pProp = iCo->pGetProperties();
        Condition::Pointer pCond =
            rReferenceCondition.Create(iCo->Id(), iCo->GetGeometry(), pProp);
        rDestinationModelPart.Conditions().push_back(pCond);
    }

    // Create new communicator local condition list
    ModelPart::ConditionsContainerType& rDestinationLocalConditions =
        pDestinationComm->LocalMesh().Conditions();
    rDestinationLocalConditions.clear();
    rDestinationLocalConditions.reserve(rDestinationModelPart.NumberOfConditions());
    for (ModelPart::ConditionsContainerType::ptr_iterator iCo =
             rDestinationModelPart.Conditions().ptr_begin();
         iCo != rDestinationModelPart.Conditions().ptr_end(); iCo++)
    {
        rDestinationLocalConditions.push_back(*iCo);
    }

    // for (auto cond: rDestinationModelPart.Conditions())
    //     KRATOS_WATCH(cond);

    // KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::CalculateYplus(unsigned int Step)
{
    int number_of_nodes = mrModelPart.NumberOfNodes();
    const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

    const double von_karman = r_current_process_info[WALL_VON_KARMAN];
    const double beta = r_current_process_info[WALL_SMOOTHNESS_BETA];

#pragma omp parallel for
    for (int i = 0; i < number_of_nodes; ++i)
    {
        NodeType& r_node = *(mrModelPart.NodesBegin() + i);

        const array_1d<double, 3>& r_velocity =
            r_node.FastGetSolutionStepValue(VELOCITY, Step);
        const double velocity_norm = norm_2(r_velocity);

        const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);

        double& y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);

        y_plus = CalculationUtilities::CalculateYplus(
            velocity_norm, wall_distance, nu, von_karman, beta, 10);
    }
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::FindConditionsParentElements(
    ModelPart* pModelPart)
{
    FindNodalNeighboursProcess find_nodal_neighbours_process(*pModelPart);
    find_nodal_neighbours_process.Execute();

    const int number_of_conditions = pModelPart->NumberOfConditions();

#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        Condition& r_condition = *(pModelPart->ConditionsBegin() + i_cond);

        Geometry<NodeType>& r_geometry = r_condition.GetGeometry();
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
            Geometry<NodeType>& r_element_geometry =
                element_candidates[i_elem].GetGeometry();
            const IndexType number_of_element_nodes = r_element_geometry.PointsNumber();
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

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::InitializeConditionFlagsForModelPart(
    ModelPart* pModelPart, const Flags& rFlag)
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

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::InitializeConditionsForModelPart(
    ModelPart* pModelPart)
{
    this->InitializeConditionFlagsForModelPart(pModelPart, INLET);
    this->InitializeConditionFlagsForModelPart(pModelPart, OUTLET);
    this->InitializeConditionFlagsForModelPart(pModelPart, STRUCTURE);

    this->FindConditionsParentElements(pModelPart);
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::AddSolutionStepVariables()
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

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::AddDofs()
{
    KRATOS_INFO("TurbulenceModel")
        << "Added eddy viscosity turbulence model dofs.\n";
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::InitializeTurbulenceModelPart()
{
    KRATOS_THROW_ERROR(std::runtime_error,
                       "Calling base class InitializeTurbulenceModelPart.", "");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::UpdateFluidViscosity()
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

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::InitializeConditions()
{
    KRATOS_THROW_ERROR(std::runtime_error,
                       "Calling base class InitializeConditions.", "");
}

/* Protected functions ****************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::InitializeNodeFlags(
    const Parameters& rParameters, const Flags& rFlag, const bool FlagValue)
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
            it_node->Set(rFlag, FlagValue);
        }
    }

    KRATOS_CATCH("");
}

/* External functions *****************************************************/

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
typedef LinearSolver<SparseSpaceType, DenseSpaceType> LinearSolverType;

template class TurbulenceEddyViscosityModelProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType>;
template class TurbulenceEddyViscosityModelProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType>;
} // namespace Kratos
