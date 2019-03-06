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
        KRATOS_INFO("TurbulenceModel") << "Finished calculating wall distances." << std::endl;

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
