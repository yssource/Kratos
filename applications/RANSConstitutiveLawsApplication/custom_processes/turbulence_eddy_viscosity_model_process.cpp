#include "turbulence_eddy_viscosity_model_process.h"

#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
/* Public functions *******************************************************/
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::TurbulenceEddyViscosityModelProcess(
    ModelPart& rModelPart, Parameters& rParameters)
    : mrModelPart(rModelPart), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "inlet_conditions"      : ["PLEASE_SPECIFY_INLET_CONDITIONS"],
        "outlet_conditions"     : ["PLEASE_SPECIFY_OUTLET_CONDITIONS"],
        "wall_conditions"       : ["PLEASE_SPECIFY_WALL_CONDITIONS"],
        "mesh_moving"       : false,
        "echo_level"        : 0,
        "model_properties"  : {},
        "distance_calculation"  : {
                "max_iterations"         : 5,
                "linear_solver_settings" : {}
        }
    })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mIsMeshMoving = mrParameters["mesh_moving"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();

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

/* Protected functions ****************************************************/
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

        y_plus = RansCalculationUtilities(this->mEchoLevel).CalculateYplus(
            velocity_norm, wall_distance, nu, von_karman, beta, 10);
    }
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


/* Protected functions ****************************************************/


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
