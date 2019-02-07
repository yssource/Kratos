#include "turbulence_eddy_viscosity_model_process.h"

namespace Kratos
{
/* Public functions *******************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::TurbulenceEddyViscosityModelProcess(
    ModelPart& rModelPart, Parameters& rParameters, TLinearSolver& rLinearSolver)
    : Process(), mrModelPart(rModelPart), mrParameters(rParameters), mrLinearSolver(rLinearSolver)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "turbulence_model_name" : "PLEASE_SPECIFY_EDDY_VISCOSITY_TURBULENCE_MODEL",
        "inlet_conditions"      : ["PLEASE_SPECIFY_INLET_CONDITIONS"],
        "outlet_conditions"     : ["PLEASE_SPECIFY_OUTLET_CONDITIONS"],
        "wall_conditions"       : ["PLEASE_SPECIFY_WALL_CONDITIONS"],
        "max_distance_calculation_iterations" : 2,
        "mesh_moving"       : False,
        "model_properties"  : {}
    })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mDistanceCalculator =
        VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
            mrModelPart, mrLinearSolver,
            mrParameters["max_distance_calculation_iterations"].GetInt());

    mIsMeshMoving = mrParameters["mesh_moving"].GetBool();

    // Add wall conditions to wall conditions list
    for (std::string model_part_name : mrParameters["wall_conditions"].GetVector())
    {
        KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(model_part_name))
            << "TurbulenceEddyViscosityModelProcess: Wall condition "
            << model_part_name << " not found." << std::endl;
        mWallConditionsList.push_back(mrModelPart.GetSubModelPart(model_part_name));
    }

    // Add inlet conditions to inlet conditions list
    for (std::string model_part_name : mrParameters["inlet_conditions"].GetVector())
    {
        KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(model_part_name))
            << "TurbulenceEddyViscosityModelProcess: Inlet condition "
            << model_part_name << " not found." << std::endl;
        mInletConditionsList.push_back(mrModelPart.GetSubModelPart(model_part_name));
    }

    // Add outlet conditions to outlet conditions list
    for (std::string model_part_name : mrParameters["outlet_conditions"].GetVector())
    {
        KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(model_part_name))
            << "TurbulenceEddyViscosityModelProcess: Outlet condition "
            << model_part_name << " not found." << std::endl;
        mOutletConditionsList.push_back(mrModelPart.GetSubModelPart(model_part_name));
    }

    // TODO: Add if blocks to navigate to correct turbulence model, this is a dummy turbulence model
    mpEddyViscosityModel = new TurbulenceEddyViscosityModel(
        mrModelPart, mrParameters, mWallConditionsList, mInletConditionsList,
        mOutletConditionsList);

    // Calculate the distances only once if the mesh is not moving
    if (!mIsMeshMoving)
        CalculateWallDistances();

    KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitialize()
{
    this->InitializeTurbulenceModelPart();

    NodesArrayType& nodes = mrTurbulenceModelPart.Nodes();

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
    {
        auto it_node = nodes.begin() + i;
        const double kinematic_viscosity = it_node->FastGetSolutionStepValue(VISCOSITY);
        it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = kinematic_viscosity;
    }

    mpEddyViscosityModel->ExecuteInitialize();
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    // Calculate the wall distances if the mesh has moved
    if (mIsMeshMoving)
        CalculateWallDistances();

    mpEddyViscosityModel->ExecuteInitializeSolutionStep();

    KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Execute()
{
    KRATOS_TRY

    mpEddyViscosityModel->Execute();

    NodesArrayType& nodes = mrTurbulenceModelPart.Nodes();

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

    // Fixing the wall boundaries for wall distance calculation
    for (ModelPart& model_part : mWallConditionsList)
    {
        NodesArrayType& nodes_array = model_part.Nodes();

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
        {
            auto it_node = nodes_array.begin() + i;
            it_node->FastGetSolutionStepValue(DISTANCE) = 0.0;
        }
    }

    mDistanceCalculator.Execute();

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
    KRATOS_TRY

    // Copy general ModelPart properites
    rDestinationModelPart.GetNodalSolutionStepVariablesList() =
        rOriginModelPart.GetNodalSolutionStepVariablesList();
    rDestinationModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());
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

    KRATOS_CATCH("");
}

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::GenerateSolutionStrategies()
{
    KRATOS_TRY

    Parameters &rParameters = *mpParameters;
    double KTolerance = rParameters["k non-linear tolerance"].GetDouble();
    double OmegaTolerance = rParameters["omega non-linear tolerance"].GetDouble();
    int KMaxIter = rParameters["maximum k non-linear iterations"].GetInt();
    int OmegaMaxIter = rParameters["maximum omega non-linear iterations"].GetInt();

    bool CalculateReactions = false;
    bool ReformDofSet = false;
    bool MoveMeshFlag = false;
    const double NearlyZero = 1.0e-8;

    // K solution strategy
    BuilderAndSolverPointerType pKBuilderAndSolver = BuilderAndSolverPointerType(
            new ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(mpKLinearSolver) );
    mpKBuilderAndSolver = pKBuilderAndSolver;

    SchemePointerType pKScheme = SchemePointerType(
            new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );

    ConvergenceCriteriaPointerType pKConvergenceCriteria = ConvergenceCriteriaPointerType(
            new DisplacementCriteria<TSparseSpace,TDenseSpace>(KTolerance,NearlyZero) );

    mpKStrategy = typename StrategyType::Pointer(
            new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            //new LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                *mpTurbulenceModelPart,
                pKScheme,
                mpKLinearSolver,
                pKConvergenceCriteria,
                pKBuilderAndSolver,
                KMaxIter,
                CalculateReactions,
                ReformDofSet,
                MoveMeshFlag)
            );

    // Omega solution strategy
    BuilderAndSolverPointerType pOmegaBuilderAndSolver = BuilderAndSolverPointerType(
            new ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver>(mpOmegaLinearSolver) );

    SchemePointerType pOmegaScheme = SchemePointerType(
            new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );

    ConvergenceCriteriaPointerType pOmegaConvergenceCriteria = ConvergenceCriteriaPointerType(
            new DisplacementCriteria<TSparseSpace,TDenseSpace>(OmegaTolerance,NearlyZero) );

    mpOmegaStrategy = typename StrategyType::Pointer(
            new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            //new LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                *mpTurbulenceModelPart,
                pOmegaScheme,
                mpOmegaLinearSolver,
                pOmegaConvergenceCriteria,
                pOmegaBuilderAndSolver,
                OmegaMaxIter,
                CalculateReactions,
                ReformDofSet,
                MoveMeshFlag)
            );

    // Convergence control for k-omega coupling iterations
    mMaximumCouplingIterations = rParameters["maximum k-omega coupling iterations"].GetInt();
    double RelativeViscosityTolerance = rParameters["turbulent viscosity tolerance"].GetDouble();
    const double AbsoluteViscosityTolerance = 1e-12; // We only want this to converge if turbulent viscosity is exactly 0 everywhere

    mpConvergenceCriteria = ConvergenceCriteriaPointerType(new TurbulentViscosityCriteria<TSparseSpace,TDenseSpace>(RelativeViscosityTolerance,AbsoluteViscosityTolerance) );
    mpConvergenceCriteria->Initialize(*mpTurbulenceModelPart);

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

} // namespace Kratos
