from __future__ import print_function, absolute_import, division

def CreateSolver(cosim_solver_settings):
    return CoSimulationBaseSolver(cosim_solver_settings)

class CoSimulationBaseSolver(object):
    """The base class for the CoSimulation Solvers
    The intention is that every solver that derives from this class
    can be used standalone.
    """
    def __init__(self, cosim_solver_settings):
        """Constructor of the Base-Solver
        Deriving classes should call it in their constructors
        """
        self.cosim_solver_settings = cosim_solver_settings

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def AdvanceInTime(self, current_time):
        pass

    def Predict(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def OutputSolutionStep(self):
        pass

    def SolveSolutionStep(self):
        pass

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