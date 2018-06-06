from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


def CreateIO(model, solver_settings):

    solver_type = solver_settings["solver_type"].GetString()

    if (solver_type == "kratos_fluid"):
        solver_module_name = "kratos_fluid_solver"
    elif (solver_type == "kratos_structural"):
        solver_module_name = "kratos_structural_solver"
    elif (solver_type == "kratos_empire"):
        solver_module_name = "kratos_empire_solver"
    else:
        raise NameError("Requested solver not available")

    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, solver_settings["solver_settings"])

    return solver
