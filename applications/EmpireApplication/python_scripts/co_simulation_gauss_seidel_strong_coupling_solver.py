from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# other imports
import co_sim_convergence_accelerator_factory as convergence_accelerator_factory
import co_sim_convergence_criteria_factory as convergence_criteria_factory
import io_factory

def CreateSolver(cosim_solver_settings):
    return GaussSeidelStrongCouplingSolver(cosim_solver_settings)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings):
        super(GaussSeidelStrongCouplingSolver, self).__init__(cosim_solver_settings)

        if not len(self.cosim_solver_settings["solvers"]) == 2:
            raise Exception("Exactly two solvers have to be specified for the GaussSeidelStrongCouplingSolver!")

        settings_solver_1 = self.cosim_solver_settings["solvers"][0]
        settings_solver_2 = self.cosim_solver_settings["solvers"][1]

        import python_solvers_wrapper_co_simulation as solvers_wrapper
        self.solver_1 = solvers_wrapper.CreateSolver(settings_solver_1)
        self.solver_2 = solvers_wrapper.CreateSolver(settings_solver_2)

    def Initialize(self):
        self.solver_1.Initialize()
        self.solver_2.Initialize()

        # Initialize FSI
        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"]
        self.convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(
            self.cosim_solver_settings["convergence_accelerator_settings"])
        self.convergence_criteria = convergence_criteria_factory.CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria_settings"])

    def Finalize(self):
        self.solver_1.FinalizeSolutionStep()
        self.solver_2.FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        new_time_1 = self.solver_1.AdvanceInTime(current_time)
        new_time_2 = self.solver_2.AdvanceInTime(current_time)
        self.convergence_accelerator.AdvanceTimeStep()

        if abs(new_time_1- new_time_2) > 1e-12:
            raise Exception("Solver time mismatch")

        return new_time_1

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
        for k in range(self.num_coupling_iterations):
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