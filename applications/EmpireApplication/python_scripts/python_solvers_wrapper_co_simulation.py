from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

available_solvers = {
    "kratos_fluid"                   : "kratos_fluid_solver",
    "kratos_structural"              : "kratos_structural_solver",
    "kratos_empire"                  : "kratos_empire_solver",
    "gauss_seidel_strong_coupling"   : "co_simulation_gauss_seidel_strong_coupling_solver",
    "weak_coupling"                  : "co_simulation_weak_coupling_solver"
}

def CreateSolver(cosim_solver_settings):

    solver_type = cosim_solver_settings["solver_type"]

    if solver_type in available_solvers:
        solver_module = __import__(available_solvers[solver_type])
        return solver_module.CreateSolver(cosim_solver_settings)
    else:
        err_msg  = 'The requested solver "' + solver_type + '" is not available!\n'
        err_msg += 'The following solvers are available:\n'
        for avail_solver in available_solvers:
            err_msg += avail_solver + "\n"
        raise NameError(err_msg)
