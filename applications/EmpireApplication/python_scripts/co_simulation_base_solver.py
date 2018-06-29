from __future__ import print_function, absolute_import, division

def CreateSolver(cosim_solver_settings):
    return CoSimulationBaseSolver(cosim_solver_settings)

class CoSimulationBaseSolver(object):
    """The base class for the CoSimulation Solvers
    The intention is that every solver that derives from this class
    can be used standalone.
    """
    def __init__(self, cosim_solver_settings):
        pass

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def AdvanceInTime(self):
        pass

    def Predict(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
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