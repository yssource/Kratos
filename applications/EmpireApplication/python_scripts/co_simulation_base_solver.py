from __future__ import print_function, absolute_import, division

def CreateSolver(model, custom_settings):
    return CoSimulationBaseSolver(model, custom_settings)

class CoSimulationBaseSolver(object):
    """The base class for the CoSimulation Solvers
    """

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