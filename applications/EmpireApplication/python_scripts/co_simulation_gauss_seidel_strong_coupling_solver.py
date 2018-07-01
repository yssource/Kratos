from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# other imports
from convergence_accelerator_factory import CreateConvergenceAccelerator
from convergence_criteria_factory import CreateConvergenceCriteria
import io_factory

def CreateSolver(cosim_solver_settings):
    return GaussSeidelStrongCouplingSolver(cosim_solver_settings)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings):
        super(GaussSeidelStrongCouplingSolver, self).__init__(cosim_solver_settings)

        settings_solver_1 = self.cosim_solver_settings["solvers"][0]
        settings_solver_2 = self.cosim_solver_settings["solvers"][1]

        import python_solvers_wrapper_co_simulation as solvers_wrapper
        self.solver_1 = solvers_wrapper.CreateSolver(settings_solver_1)
        self.solver_2 = solvers_wrapper.CreateSolver(settings_solver_2)


    def Initialize(self):
        self.solver_1.Initialize()
        self.solver_2.Initialize()

        # Initialize FSI
        self.num_coupling_iterations = self.settings["num_coupling_iterations"]
        self.convergence_accelerator = CreateConvergenceAccelerator(self.settings["convergence_accelerator_settings"])
        self.convergence_criteria = CreateConvergenceCriteria(self.settings["convergence_criteria_settings"])
        # TODO check if the time settings are consistent!

    def Finalize(self):
        self.solver_1.FinalizeSolutionStep()
        self.solver_2.FinalizeSolutionStep()

    def AdvanceInTime(self):
        self.solver_1.AdvanceInTime()
        self.solver_2.AdvanceInTime()
        self.convergence_accelerator.AdvanceTimeStep()
        # TODO return time and check if it is consistent!

    def Predict(self):
        self.solver_1.Predict()
        self.solver_2.Predict()

    def InitializeSolutionStep(self):
        self.solver_1.InitializeSolutionStep()
        self.solver_2.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.solver_1.FinalizeSolutionStep()
        self.solver_2.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for k in range(max_iter):
            self.solver_1.ImportData()
            self.solver_1.SolveSolutionStep()
            self.solver_1.ExportData()
            self.solver_2.ImportData()
            self.solver_2.SolveSolutionStep()
            self.solver_2.ExportData()
            if self.convergence_criteria.IsConverged():
                break
            else:
                self.convergence_accelerator.ComputeUpdate(...)


    def ImportData(self, DataName, FromClient):
        raise NotImplementedError("This needs to be implemented!")
    def ImportMesh(self, MeshName, FromClient):
        raise NotImplementedError("This needs to be implemented!")

    def ExportData(self, DataName, ToClient):
        raise NotImplementedError("This needs to be implemented!")
    def ExportMesh(self, MeshName, ToClient):
        raise NotImplementedError("This needs to be implemented!")

    def MakeDataAvailable(self, DataName, ToClient):
        raise NotImplementedError("This needs to be implemented!")
    def MakeMeshAvailable(self, MeshName, ToClient):
        raise NotImplementedError("This needs to be implemented!")