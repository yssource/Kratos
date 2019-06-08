from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# imports Kratos tools
import co_simulation_tools as cs_tools

# import Kratos components
from co_simulation_predictors.co_simulation_predictor_factory import CreatePredictor
from co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory import CreateConvergenceAccelerator
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria

# import Kratos solvers
from co_simulation_solvers.sdof_solver import SDoFSolver
from co_simulation_solvers.kratos_fluid_solver import KratosFluidSolver

### WORK HERE: imports to be resolved
import co_simulation_solvers.python_solvers_wrapper_co_simulation as solvers_wrapper

# imports from basic Python
import json
import sys




# read parameter file(s)
parameter_file_name = "parameters_FSI_AllInOne.json"
with open(parameter_file_name, 'r') as parameter_file:
    cosim_settings = json.load(parameter_file)


### CoSimulationAnalysis - init ###################################
if type(cosim_settings) != dict:
    raise Exception("Input is expected to be provided as a python dictionary")

# commented out to enable changes in the settings structure without changing cs_tools
# cs_tools.CheckCoSimulationSettingsAndAssignDefaults(cosim_settings)

problem_data = cosim_settings["problem_data"]

cs_tools.PRINT_COLORS = problem_data["print_colors"]

parallel_type = problem_data["parallel_type"]
if parallel_type == "OpenMP":
    flush_stdout = True
    cs_tools.PRINTING_RANK = True
elif parallel_type == "MPI":
    from co_simulation_mpi_space import CoSimulationMPISpace

    cs_tools.COSIM_SPACE = CoSimulationMPISpace()
    cs_tools.PRINTING_RANK = (cs_tools.COSIM_SPACE.Rank() == 0)
else:
    raise Exception('"parallel_type" can only be "OpenMP" or "MPI"!')

flush_stdout = problem_data["flush_stdout"]
echo_level = problem_data["echo_level"]


### Initialize ###################################
cosim_solver_settings = cosim_settings["solver_settings"]

## CoSimulationBaseSolver __init__() ###################################
echo_level = 0
solver_level = 0
if "echo_level" in cosim_solver_settings:
    echo_level = cosim_solver_settings["echo_level"]
io_is_initialized = False

## CoSimulationBaseCouplingSolver __init__() ###################################
solver_names = ["fluid", "structure"]
solvers = {}

solvers["fluid"] = KratosFluidSolver(cosim_solver_settings["solvers"]["fluid"], solver_level-1)
solvers["structure"] = SDoFSolver(cosim_solver_settings["solvers"]["structure"], solver_level-1)

if "predictor_settings" in cosim_solver_settings:
    predictor = CreatePredictor(cosim_solver_settings["predictor_settings"], solvers, solver_level)
    predictor.SetEchoLevel(echo_level)
else:
    predictor = None

# With this setting the coupling can start later
start_coupling_time = 0.0
if "start_coupling_time" in cosim_solver_settings:
    start_coupling_time = cosim_solver_settings["start_coupling_time"]
if start_coupling_time > 0.0:
    coupling_started = False
else:
    coupling_started = True

## GaussSeidelStrongCouplingSolver __init__() ###################################
convergence_accelerator = CreateConvergenceAccelerator(
    cosim_solver_settings["convergence_accelerator_settings"], solvers, solver_level)
convergence_accelerator.SetEchoLevel(echo_level)

convergence_criteria = CreateConvergenceCriteria(
    cosim_solver_settings["convergence_criteria_settings"], solvers, solver_level)
convergence_criteria.SetEchoLevel(echo_level)

num_coupling_iterations = cosim_solver_settings["num_coupling_iterations"]

## GaussSeidelStrongCouplingSolver Initialize() ###################################
# from super -------------------
for solver_name in solver_names:
    solvers[solver_name].Initialize()
for solver_name in solver_names:
    solvers[solver_name].InitializeIO(solvers, echo_level)

if predictor is not None:
    predictor.Initialize()
# from sub -------------------
convergence_accelerator.Initialize()
convergence_criteria.Initialize()

## GaussSeidelStrongCouplingSolver Check() ###################################
# from super -------------------
for solver_name in solver_names:
    solvers[solver_name].Check()

if predictor is not None:
    predictor.Check()

# from sub -------------------
convergence_accelerator.Check()
convergence_criteria.Check()


if echo_level > 0:
    ## GaussSeidelStrongCouplingSolver PrintInfo() ###################################
    # from super -------------------
    cs_tools.cs_tools.couplingsolverprint("GaussSeidelStrongCouplingSolver", "Has the following participants:")
    for solver_name in solver_names:
        solvers[solver_name].PrintInfo()

    if predictor is not None:
        cs_tools.cs_tools.couplingsolverprint("GaussSeidelStrongCouplingSolver", "Uses a Predictor:")
        predictor.PrintInfo()

    # from sub -------------------
    cs_tools.cs_tools.couplingsolverprint("GaussSeidelStrongCouplingSolver", "Uses the following objects:")
    convergence_accelerator.PrintInfo()
    convergence_criteria.PrintInfo()

# Stepping and time settings
end_time = cosim_settings["problem_data"]["end_time"]
time = cosim_settings["problem_data"]["start_time"]
step = 0

if flush_stdout:
    sys.stdout.flush()

### RunSolutionLoop ###################################
print("")
while time < end_time:
    print("")
    step += 1

    ## GaussSeidelStrongCouplingSolver AdvanceInTime() ###################################
    # from super -------------------
    old_time = time
    time = solvers[solver_names[0]].AdvanceInTime(old_time)
    for solver_name in solver_names[1:]:
        time_other_solver = solvers[solver_name].AdvanceInTime(old_time)
        if abs(time - time_other_solver) > 1e-12:
            raise Exception("Solver time mismatch")

    if not coupling_started and time > start_coupling_time:
        coupling_started = True
        if echo_level > 0:
            cs_tools.couplingsolverprint(solver_level, "GaussSeidelStrongCouplingSolver", cs_tools.bold("Starting Coupling"))

    # if a predictor is used then the delta_time is set
    # this is needed by some predictors
    if predictor is not None:
        delta_time = time - old_time
        predictor.SetDeltaTime(delta_time)

    cs_tools.csprint(0, cs_tools.bold("time={0:.12g}".format(time) + " | step=" + str(step)))
    ## GaussSeidelStrongCouplingSolver InitializeSolutionStep() ###################################
    # from super -------------------
    for solver_name in solver_names:
        solvers[solver_name].InitializeSolutionStep()

    if predictor is not None:
        predictor.InitializeSolutionStep()
    convergence_accelerator.InitializeSolutionStep()
    convergence_criteria.InitializeSolutionStep()

    ## GaussSeidelStrongCouplingSolver SolveSolutionStep() ###################################
    for k in range(num_coupling_iterations):
        if echo_level > 0:
            cs_tools.couplingsolverprint(solver_level, "GaussSeidelStrongCouplingSolver",
                                cs_tools.cyan("Coupling iteration:"),
                                cs_tools.bold(str(k + 1) + " / " + str(num_coupling_iterations)))

        convergence_accelerator.InitializeNonLinearIteration()
        convergence_criteria.InitializeNonLinearIteration()

        for solver_name in solver_names:
            solver = solvers[solver_name]

            ## GaussSeidelStrongCouplingSolver _SynchronizeInputData() ###################################
            if coupling_started:
                input_data_list = cosim_solver_settings["coupling_loop"][solver_name]["input_data_list"]

                if time >= cosim_solver_settings["coupling_loop"][solver_name]["input_coupling_start_time"]:
                    for input_data in input_data_list:
                        from_solver = solvers[input_data["from_solver"]]
                        data_name = input_data["data_name"]
                        data_definition = from_solver.GetDataDefinition(data_name)
                        data_settings = {"data_format": data_definition["data_format"],
                                         "data_name": data_name,
                                         "io_settings": input_data["io_settings"]}
                        solver.ImportData(data_settings, from_solver)

            solver.SolveSolutionStep()

            ## GaussSeidelStrongCouplingSolver _SynchronizeOutputData() ###################################
            if coupling_started:
                output_data_list = cosim_solver_settings["coupling_loop"][solver_name]["output_data_list"]

                if time >= cosim_solver_settings["coupling_loop"][solver_name]["output_coupling_start_time"]:
                    for output_data in output_data_list:
                        to_solver = solvers[output_data["to_solver"]]
                        data_name = output_data["data_name"]
                        data_definition = to_solver.GetDataDefinition(data_name)
                        data_settings = {"data_format": data_definition["data_format"],
                                         "data_name": data_name,
                                         "io_settings": output_data["io_settings"]}
                        solver.ExportData(data_settings, to_solver)

        convergence_accelerator.FinalizeNonLinearIteration()
        convergence_criteria.FinalizeNonLinearIteration()

        if convergence_criteria.IsConverged():
            if echo_level > 0:
                cs_tools.couplingsolverprint(solver_level, "GaussSeidelStrongCouplingSolver", cs_tools.green("### CONVERGENCE WAS ACHIEVED ###"))
            break
        else:
            convergence_accelerator.ComputeUpdate()

        if k + 1 >= num_coupling_iterations and echo_level > 0:
            cs_tools.couplingsolverprint(solver_level, "GaussSeidelStrongCouplingSolver", cs_tools.red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))

    ## GaussSeidelStrongCouplingSolver FinalizeSolutionStep() ###################################
    # from super -------------------
    for solver_name in solver_names:
        solvers[solver_name].FinalizeSolutionStep()
        if predictor is not None:
            predictor.FinalizeSolutionStep()

    # from sub -------------------
    convergence_accelerator.FinalizeSolutionStep()
    convergence_criteria.FinalizeSolutionStep()


    ## GaussSeidelStrongCouplingSolver OutputSolutionStep() ###################################
    for solver_name in solver_names:
        solvers[solver_name].OutputSolutionStep()
    
    if flush_stdout:
        sys.stdout.flush()

## GaussSeidelStrongCouplingSolver Finalize() ###################################
for solver_name in solver_names:
    solvers[solver_name].Finalize()

if predictor is not None:
    predictor.Finalize()



