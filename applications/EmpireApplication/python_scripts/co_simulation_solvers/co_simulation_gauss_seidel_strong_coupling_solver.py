from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_solvers.co_simulation_base_coupling_solver import CoSimulationBaseCouplingSolver

# Other imports
from co_simulation_convergence_accelerators.co_simulation_convergence_accelerator_factory import CreateConvergenceAccelerator
from co_simulation_convergence_criteria.co_simulation_convergence_criteria_factory import CreateConvergenceCriteria
from co_simulation_tools import csprint, red, green, cyan, bold

def CreateSolver(cosim_solver_settings, level):
    return GaussSeidelStrongCouplingSolver(cosim_solver_settings, level)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, cosim_solver_settings, level):
        if not len(cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(GaussSeidelStrongCouplingSolver, self).__init__(cosim_solver_settings, level)

    def Initialize(self):
        super(GaussSeidelStrongCouplingSolver, self).Initialize()

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]

        self.convergence_accelerator = CreateConvergenceAccelerator(
            self.cosim_solver_settings["convergence_accelerator_settings"],
            self.solvers, self.cosim_solver_details, self.lvl)

        self.convergence_criteria = CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"],
            self.solvers, self.cosim_solver_details, self.lvl)

    def FinalizeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).FinalizeSolutionStep()
        self.convergence_accelerator.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            csprint(self.lvl, cyan("Coupling iteration: ")+bold(str(k+1)+" / " + str(self.num_coupling_iterations)))
            self.convergence_accelerator.SetPreviousSolution()
            self.convergence_criteria.SetPreviousSolution()
            for solver_name in self.solver_names:
                solver = self.solvers[solver_name]
                self._SynchronizeInputData(solver, solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver, solver_name)

            if self.convergence_criteria.IsConverged():
                csprint(self.lvl, green("##### CONVERGENCE AT INTERFACE WAS ACHIEVED #####"))
                break
            else:
                self.convergence_accelerator.ComputeUpdate()

            if k+1 >= self.num_coupling_iterations:
                csprint(self.lvl, red("XXXXX CONVERGENCE AT INTERFACE WAS NOT ACHIEVED XXXXX"))