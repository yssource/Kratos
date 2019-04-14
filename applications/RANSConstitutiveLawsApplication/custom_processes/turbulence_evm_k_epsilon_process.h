//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_TURBULENCE_EVM_K_EPSILON_PROCESS_H_INCLUDED)
#define KRATOS_TURBULENCE_EVM_K_EPSILON_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// Project includes
#include "includes/cfd_variables.h"
#include "includes/communicator.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/color_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_processes/turbulence_eddy_viscosity_model_process.h"
#include "custom_strategies/evm_k_epsilon/residual_based_bossak_turbulent_energy_dissipation_scheme.h"
#include "custom_strategies/evm_k_epsilon/residual_based_bossak_turbulent_kinetic_energy_scheme.h"
#include "custom_strategies/general_convergence_criteria.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_constitutive_laws_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class TurbulenceEvmKEpsilonProcess
    : public TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TurbulenceEvmKEpsilonProcess
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEvmKEpsilonProcess);

    typedef TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> SolvingStrategyType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> StrategyType;

    typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

    typedef typename SchemeType::Pointer SchemePointerType;

    typedef typename SchemeType::DofsArrayType DofsArrayType;

    typedef typename SchemeType::TSystemMatrixType SystemMatrixType;

    typedef typename SchemeType::TSystemVectorType SystemVectorType;

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BuilderAndSolverType;

    typedef typename BuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> ConvergenceCriteriaType;

    typedef typename ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TurbulenceEvmKEpsilonProcess(ModelPart& rModelPart,
                                 Parameters& rParameters,
                                 typename TLinearSolver::Pointer pKLinearSolver,
                                 typename TLinearSolver::Pointer pEpsilonLinearSolver)
        : TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
              rModelPart, rParameters),
          mpKLinearSolver(pKLinearSolver),
          mpEpsilonLinearSolver(pEpsilonLinearSolver)
    {
        KRATOS_TRY

        KRATOS_INFO("TurbulenceModel")
            << "Initializing k-epsilon turbulence model" << std::endl;

        Parameters default_parameters = Parameters(R"(
        {
            "turbulent_kinetic_energy_settings": {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "max_iterations": 10,
                "echo_level": 0,
                "linear_solver_settings": {}
            },
            "turbulent_energy_dissipation_rate_settings": {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "max_iterations": 10,
                "echo_level": 0,
                "linear_solver_settings": {}
            },
            "convergence_tolerances":
            {
                "turbulent_viscosity_relative_tolerance": 1e-3,
                "turbulent_viscosity_absolute_tolerance": 1e-5,
                "maximum_coupling_iterations": 10,
                "maximum_stabilization_multiplier":1e+4,
                "echo_level": 2
            },
            "echo_level"     : 0,
            "time_scheme"    : "steady",
            "scheme_settings": {},
            "constants":
            {
                "wall_smoothness_beta"    : 5.2,
                "von_karman"              : 0.41,
                "c_mu"                    : 0.09,
                "c1"                      : 1.44,
                "c2"                      : 1.92,
                "sigma_k"                 : 1.0,
                "sigma_epsilon"           : 1.3
            },
            "flow_parameters":
            {
                "ramp_up_time"                : 0.03,
                "turbulent_viscosity_min"     : 1e-6,
                "turbulent_viscosity_max"     : 1e+2
            }
        })");

        this->mrParameters["model_properties"].ValidateAndAssignDefaults(default_parameters);
        this->mrParameters["model_properties"]["convergence_tolerances"].ValidateAndAssignDefaults(
            default_parameters["convergence_tolerances"]);
        this->mrParameters["model_properties"]["constants"].ValidateAndAssignDefaults(
            default_parameters["constants"]);
        this->mrParameters["model_properties"]["flow_parameters"].ValidateAndAssignDefaults(
            default_parameters["flow_parameters"]);

        const Parameters& r_convergence_parameters =
            this->mrParameters["model_properties"]["convergence_tolerances"];

        this->mTurbulentViscosityRelativeTolerance =
            r_convergence_parameters["turbulent_viscosity_relative_tolerance"].GetDouble();
        this->mTurbulentViscosityAbsoluteTolerance =
            r_convergence_parameters["turbulent_viscosity_absolute_tolerance"].GetDouble();
        this->mMaximumCouplingIterations =
            r_convergence_parameters["maximum_coupling_iterations"].GetInt();


        if (this->mrParameters["model_properties"]["time_scheme"].GetString() ==
            "steady")
        {
            rModelPart.GetProcessInfo()[RANS_TIME_STEP] = 1;
        }
        else if (this->mrParameters["model_properties"]["time_scheme"].GetString() ==
                 "transient")
        {
            rModelPart.GetProcessInfo()[RANS_TIME_STEP] = 0;
        }
        else
        {
            KRATOS_ERROR
                << "Undefined time_scheme in turbulence_model_settings. Only "
                   "\"steady\" or \"transient\" is allowed in k_epsilon model "
                   "time_scheme.";
        }

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~TurbulenceEvmKEpsilonProcess()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetIsCoSolvingProcessActive(bool rIsSolvingProcessActive)
    {
        this->IsSolvingProcessActive = rIsSolvingProcessActive;
    }

    void AddStrategy(SolvingStrategyType* pStrategy)
    {
        mrSolvingStrategiesList.push_back(pStrategy);
    }

    virtual void Execute() override
    {
        if (!this->IsSolvingProcessActive)
            return;

        KRATOS_INFO_IF("TurbulenceModel", this->mEchoLevel > 0)
            << "Solving RANS equations...\n";
        this->SolveStep();
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        BaseType::ExecuteInitializeSolutionStep();

        if (!this->IsSolvingProcessActive)
            return;

        AssignInitialConditions();
        this->UpdateBoundaryConditions();

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY

        BaseType::ExecuteFinalizeSolutionStep();

        this->UpdateFluidViscosity();

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
    virtual std::string Info() const override
    {
        return std::string("TurbulenceEvmKEpsilonProcess");
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::map<std::size_t, double> mRatioTolerance;
    std::map<std::size_t, double> mAbsTolerance;

    typename TLinearSolver::Pointer mpKLinearSolver;

    typename StrategyType::Pointer mpKStrategy;

    SchemePointerType mpKScheme;

    SchemePointerType mpEpsilonScheme;

    typename TLinearSolver::Pointer mpEpsilonLinearSolver;

    typename StrategyType::Pointer mpEpsilonStrategy;

    double mFreestreamK;

    int mMaximumCouplingIterations;

    double mTurbulentViscosityRelativeTolerance;

    double mTurbulentViscosityAbsoluteTolerance;

    bool mAssignedInitialConditions = false;

    ConvergenceCriteriaPointerType mpConvergenceCriteria;

    BuilderAndSolverPointerType mpKBuilderAndSolver;

    BuilderAndSolverPointerType mpEpsilonBuilderAndSolver;

    ModelPart* mpTurbulenceKModelPart;

    ModelPart* mpTurbulenceEpsilonModelPart;

    bool IsSolvingProcessActive = false;

    std::vector<SolvingStrategyType*> mrSolvingStrategiesList;

    ///@}
    ///@name Private Operators
    ///@{

    void AssignInitialConditions()
    {
        if (mAssignedInitialConditions)
            return;

        this->CalculateYplus();
        this->UpdateInletNodes();

        const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();
        unsigned int initialized_nodes_count{0};

        const int number_of_nodes = this->mrModelPart.NumberOfNodes();
#pragma omp parallel for reduction(+ : initialized_nodes_count)
        for (int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodeIterator i_node = this->mrModelPart.NodesBegin() + i;
            if (!i_node->IsFixed(TURBULENT_KINETIC_ENERGY) &&
                !i_node->IsFixed(TURBULENT_ENERGY_DISSIPATION_RATE))
            {
                double k{0.0}, epsilon{0.0};
                this->CalculateTurbulentValues(k, epsilon, *i_node, r_current_process_info);
                this->InitializeValues(*i_node, k, epsilon);
                initialized_nodes_count++;
            }
        }

        KRATOS_INFO("TurbulenceModel")
            << initialized_nodes_count << " nodes were assigned initial values in "
            << this->mrModelPart.Name() << " model part.\n";

        mAssignedInitialConditions = true;
    }

    void InitializeValues(NodeType& rNode, double K, double Epsilon)
    {
        rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = K;
        rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = Epsilon;
    }

    void UpdateBoundaryConditions()
    {
        this->CalculateYplus();
        this->UpdateInletNodes();
        this->UpdateWallNodes();
    }

    void UpdateInletNodes()
    {
        const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

        // Set Boundary Conditions
        const int number_of_nodes = this->mrModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodeIterator i_node = this->mrModelPart.NodesBegin() + i;
            if (i_node->Is(INLET))
            {
                double k{0.0}, epsilon{0.0};
                this->CalculateTurbulentValues(k, epsilon, *i_node, r_current_process_info);
                this->InitializeValues(*i_node, k, epsilon);
            }
        }
    }

    void UpdateWallNodes()
    {
        // Set Boundary Conditions
        const int number_of_nodes = this->mrModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodeIterator i_node = this->mrModelPart.NodesBegin() + i;
            if (i_node->Is(STRUCTURE))
            {
                this->InitializeValues(*i_node, 1e-10, 1e-10);
            }
        }
    }

    void SolveStep()
    {
        this->UpdateTurbulentViscosity();

        for (auto strategy: mrSolvingStrategiesList)
        {
            strategy->InitializeSolutionStep();
            strategy->Predict();
        }

        const int number_of_nodes = this->mrModelPart.NumberOfNodes();
        Vector old_turbulent_viscosity = ZeroVector(number_of_nodes);
        Vector new_turbulent_viscosity = ZeroVector(number_of_nodes);

        bool is_converged = false;
        int step = 1;

        DofsArrayType dummy_dofs;
        SystemMatrixType dummy_matrix;
        SystemVectorType dummy_vector;

        while (!is_converged && (step <= mMaximumCouplingIterations))
        {
#pragma omp parallel for
            for (int i = 0; i < number_of_nodes; ++i)
            {
                old_turbulent_viscosity[i] =
                    (this->mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            }

            for (auto strategy: mrSolvingStrategiesList)
                strategy->SolveSolutionStep();

            this->UpdateTurbulentViscosity();

#pragma omp parallel for
            for (int i = 0; i < number_of_nodes; ++i)
            {
                new_turbulent_viscosity[i] =
                    (this->mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            }

            double increase_norm(0.0), solution_norm(0.0);

#pragma omp parallel for reduction(+ : increase_norm, solution_norm)
            for (int i = 0; i < number_of_nodes; ++i)
            {
                increase_norm += new_turbulent_viscosity[i] - old_turbulent_viscosity[i];
                solution_norm += new_turbulent_viscosity[i];
            }

            if (solution_norm == 0.0)
                solution_norm = 1.0;

            const double ratio = std::abs(increase_norm / solution_norm);
            increase_norm = std::abs(increase_norm) / number_of_nodes;

            if (this->mEchoLevel > 0)
                std::cout << "[" << step << "] CONVERGENCE CHECK: TURBULENT_VISCOSITY: ratio = "
                          << std::scientific << ratio << "; exp.ratio = " << std::scientific
                          << mTurbulentViscosityRelativeTolerance
                          << ": abs = " << std::scientific << increase_norm
                          << "; exp.abs = " << std::scientific
                          << mTurbulentViscosityAbsoluteTolerance << std::endl;

            is_converged = (increase_norm < mTurbulentViscosityAbsoluteTolerance) ||
                           (ratio < mTurbulentViscosityRelativeTolerance);

            if (this->mEchoLevel > 0 && is_converged)
                std::cout << "[" << step << "] "
                          << "CONVERGENCE CHECK: TURBULENT_VISCOSITY: "
                             "*** CONVERGENCE IS ACHIEVED ***\n";

            step++;
        }

        if (!is_converged)
        {
            std::cout
                << "|----------------------------------------------------|"
                << std::endl;
            std::cout << "|    " << BOLDFONT(FRED("ATTENTION: Max coupling iterations exceeded"))
                      << "     |" << std::endl;
            std::cout
                << "|----------------------------------------------------|"
                << std::endl;
        }

            for (auto strategy: mrSolvingStrategiesList)
                strategy->FinalizeSolutionStep();

    }

    double CalculateNodalTurbulentViscosity(const ModelPart::NodeIterator& iNode,
                                            const double C_mu,
                                            const unsigned int Step) const
    {
        const double tke = iNode->FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY, Step);
        const double epsilon =
            iNode->FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE, Step);
        const double y_plus = iNode->FastGetSolutionStepValue(RANS_Y_PLUS);
        const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);

        const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
            C_mu, tke, epsilon, f_mu);

        return nu_t;
    }

    void UpdateTurbulentViscosity()
    {
        DofsArrayType dummy_dofs;
        SystemMatrixType dummy_matrix;
        SystemVectorType dummy_vector;

        const int number_of_nodes = this->mrModelPart.NumberOfNodes();

        const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

        const double C_mu = r_current_process_info[TURBULENCE_RANS_C_MU];
        const double nu_t_min = r_current_process_info[TURBULENT_VISCOSITY_MIN];
        const double nu_t_max = r_current_process_info[TURBULENT_VISCOSITY_MAX];

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodeIterator iNode = this->mrModelPart.NodesBegin() + i;
            const double nu_t = CalculateNodalTurbulentViscosity(iNode, C_mu, 0);
            iNode->FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
        }

        RansCalculationUtilities(this->mEchoLevel)
            .WarnIfNegative(this->mrModelPart, TURBULENT_VISCOSITY, "NuT");
        RansCalculationUtilities(this->mEchoLevel)
            .ClipVariable(this->mrModelPart, TURBULENT_VISCOSITY, nu_t_min, nu_t_max);
    }

    void CalculateTurbulentValues(double& TurbulentKineticEnergy,
                                  double& TurbulentEnergyDissipationRate,
                                  const NodeType& rNode,
                                  const ProcessInfo& rProcessInfo)
    {
        const double c_mu = rProcessInfo[TURBULENCE_RANS_C_MU];
        const double von_karman = rProcessInfo[WALL_VON_KARMAN];
        const double y_plus = rNode.FastGetSolutionStepValue(RANS_Y_PLUS, 0);
        const double nu = rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double wall_distance = rNode.FastGetSolutionStepValue(DISTANCE);

        EvmKepsilonModelUtilities::CalculateTurbulentValues(
            TurbulentKineticEnergy, TurbulentEnergyDissipationRate, y_plus, nu,
            wall_distance, c_mu, von_karman);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TurbulenceEvmKEpsilonProcess& operator=(TurbulenceEvmKEpsilonProcess const& rOther);

    /// Copy constructor.
    TurbulenceEvmKEpsilonProcess(TurbulenceEvmKEpsilonProcess const& rOther);

    ///@}

}; // Class TurbulenceEvmKEpsilonProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator>>(
    std::istream& rIStream,
    TurbulenceEvmKEpsilonProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEvmKEpsilonProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_TURBULENCE_EVM_K_EPSILON_PROCESS_H_INCLUDED  defined
