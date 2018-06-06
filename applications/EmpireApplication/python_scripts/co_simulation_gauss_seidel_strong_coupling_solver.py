from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# other imports
from convergence_accelerator_factory import CreateConvergenceAccelerator
from convergence_criteria_factory import CreateConvergenceCriteria

def CreateSolver(model, custom_settings):
    return GaussSeidelStrongCouplingSolver(model, custom_settings)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseSolver):
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, settings):
        """The constructor of the PythonSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedSolver, self).__init__(settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        settings -- The solver settings used
        """
        pass

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

    def AdvanceInTime(self, current_time):
        self.solver_1.AdvanceInTime(current_time)
        self.solver_2.AdvanceInTime(current_time)
        self.convergence_accelerator.AdvaneTimeStep()
        # TODO return time

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
        pass
    def ImportMesh(self, MeshName, FromClient):
        pass

    def ExportData(self, DataName, ToClient):
        pass
    def ExportMesh(self, MeshName, ToClient):
        pass

    def MakeDataAvailable(self, DataName, ToClient):
        pass
    def MakeMeshAvailable(self, MeshName, ToClient):
        pass