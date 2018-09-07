from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class of the cosim-solvers
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver

available_solvers = {
    "dummy"                        : "co_simulation_base_solver",
    "kratos_fluid"                 : "kratos_fluid_solver",
    "kratos_structural"            : "kratos_structural_solver",
    "kratos_empire"                : "kratos_empire_solver",
    "gauss_seidel_strong_coupling" : "co_simulation_gauss_seidel_strong_coupling_solver",
    "gauss_seidel_weak_coupling"   : "co_simulation_gauss_seidel_weak_coupling_solver",
    "sdof"                         : "sdof_solver",
    "mdof"                         : "mdof_solver"
}

def CreateSolver(cosim_solver_settings, level):
    """This function creates and returns the solvers used for CoSimulation
    New solvers have to be registered by adding them to "available_solvers"
    """
    if (type(cosim_solver_settings) != dict):
        raise Exception("Input is expected to be provided as a python dictionary")

    solver_type = cosim_solver_settings["solver_type"]

    if solver_type in available_solvers:
        solver_module = __import__(available_solvers[solver_type])
        solver = solver_module.CreateSolver(cosim_solver_settings, level+1)
        if not isinstance(solver, CoSimulationBaseSolver):
            err_msg  = 'The requested solver "' + solver_type
            err_msg += '" does not derive from "CoSimulationBaseSolver"!'
            raise Exception(err_msg)
        return solver
    else:
        err_msg  = 'The requested solver "' + solver_type + '" is not available!\n'
        err_msg += 'The following solvers are available:\n'
        for avail_solver in available_solvers:
            err_msg += "\t" + avail_solver + "\n"
        raise NameError(err_msg)
