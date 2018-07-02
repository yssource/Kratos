def GetSolverCoSimulationDetails(co_simulation_solver_settings):
    num_solvers = len(co_simulation_solver_settings)
    solver_cosim_details = {}
    for solver_settings in co_simulation_solver_settings:
        solver_name = solver_settings["name"]
        solver_cosim_details[solver_name] = solver_settings

    return solver_cosim_details

