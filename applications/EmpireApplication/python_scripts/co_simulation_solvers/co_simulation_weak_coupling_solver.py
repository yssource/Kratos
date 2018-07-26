from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_solvers.co_simulation_base_coupling_solver import CoSimulationBaseCouplingSolver

# Other imports
from co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory import CreateConvergenceAccelerator
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from co_simulation_tools import csprint, red, green, cyan, bold

def CreateSolver(cosim_solver_settings, level):
    return WeakCouplingSolver(cosim_solver_settings, level)

class WeakCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, cosim_solver_settings, level):
        if not len(cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(WeakCouplingSolver, self).__init__(cosim_solver_settings, level)

    def Initialize(self):
        super(WeakCouplingSolver, self).Initialize()

    def FinalizeSolutionStep(self):
        super(WeakCouplingSolver, self).FinalizeSolutionStep()

    def SolveSolutionStep(self):

        for solver_name in self.solver_names:
            solver = self.solvers[solver_name]
            self._SynchronizeInputData(solver, solver_name)
            solver.SolveSolutionStep()
            self._SynchronizeOutputData(solver, solver_name)

        csprint(self.lvl, green("##### SOLVED SYSTEM #####"))
