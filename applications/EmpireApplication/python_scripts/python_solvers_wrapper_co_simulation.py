from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics


def CreateSolver(model, solver_settings):

    if (type(model) != KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")

    if (type(solver_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    parallelism = solver_settings["problem_data"]["parallel_type"].GetString()
    solver_type = solver_settings["solver_settings"]["solver_type"].GetString()

    if (solver_type == "co_simulation_gauss_seidel_strong_coupling_solver"):
        solver_module_name = "co_simulation_gauss_seidel_strong_coupling_solver"
    else:
        raise NameError("Requested solver not available")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, solver_settings["solver_settings"])

    return solver
