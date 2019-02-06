#include "turbulence_evm_k_epsilon_process.h"

namespace Kratos
{
/* Public functions *******************************************************/
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::TurbulenceEvmKEpsilonProcess(
    ModelPart& rModelPart, Parameters& rParameters, TLinearSolver& rLinearSolver)
    : Process(), mrModelPart(rModelPart), mrParameters(rParameters), mrLinearSolver(rLinearSolver)
{
    Parameters default_parameters = Parameters(R"(
    {
        "inlet_conditions"  : ["PLEASE_SPECIFY_INLET_BOUNDARIES"],
        "outlet_conditions" : ["PLEASE_SPECIFY_OUTLET_BOUNDARIES"],
        "wall_conditions"   : ["PLEASE_SPECIFY_WALL_BOUNDARIES"],
        "max_distance_calculation_iterations" : 2
    })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mDistanceCalculator =
        VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
            mrModelPart, mrLinearSolver,
            mrParameters["max_distance_calculation_iterations"].GetInt());
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::~TurbulenceEvmKEpsilonProcess()
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::Execute()
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitialize()
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitializeSolutionStep()
{

}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
std::string TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::Info() const
{
    return "TurbulenceEvmKEpsilonProcess";
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << "TurbulenceEvmKEpsilonProcess";
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::PrintData(
    std::ostream& rOStream) const
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>::CalculateWallDistances()
{
    KRATOS_TRY

    // Fixing the wall boundaries for wall distance calculation
    for(auto model_part_name: mrParameters["wall_conditions"].GetVector())
    {
        KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(model_part_name))<<"TurbulenceModel: "<<model_part_name<<" not found."<<std::endl;

        NodesArrayType& nodes_array = mrModelPart.GetSubModelPart(model_part_name).Nodes();

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

/* External functions *****************************************************/

/// output stream function
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos
