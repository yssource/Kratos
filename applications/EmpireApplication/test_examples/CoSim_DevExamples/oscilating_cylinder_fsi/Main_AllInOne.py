from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# imports in CoSimulationAnalysis
import co_simulation_tools as cs_tools
from co_simulation_tools import csprint, bold, CheckCoSimulationSettingsAndAssignDefaults
import co_simulation_solvers.python_solvers_wrapper_co_simulation as solvers_wrapper

# imports from basic Python
import json
import sys

# read parameter file(s)
parameter_file_name = "project_parameters_cosim_oscilating_cylinder_fsi.json"
with open(parameter_file_name, 'r') as parameter_file:
    cosim_settings = json.load(parameter_file)


### CoSimulationAnalysis - init ###################################
if type(cosim_settings) != dict:
    raise Exception("Input is expected to be provided as a python dictionary")

CheckCoSimulationSettingsAndAssignDefaults(cosim_settings)

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
solver = solvers_wrapper.CreateSolver(cosim_settings["solver_settings"], level=0)
solver.Initialize()
solver.Check()

if echo_level > 0:
    solver.PrintInfo()

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
    time = solver.AdvanceInTime(time)

    ## InitializeSolutionStep ###################################
    csprint(0, bold("time={0:.12g}".format(time) + " | step=" + str(step)))
    solver.InitializeSolutionStep()

    solver.Predict()
    solver.SolveSolutionStep()

    ## FinalizeSolutionStep ###################################
    solver.FinalizeSolutionStep()

    ## OutputSolutionStep ###################################
    solver.OutputSolutionStep()

    if flush_stdout:
        sys.stdout.flush()

    solver.Finalize()



################################### BAUSTELLE ###################################
