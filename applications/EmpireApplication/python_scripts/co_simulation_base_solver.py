from __future__ import print_function, absolute_import, division

# Other imports
import io_factory

def CreateSolver(cosim_solver_settings, level):
    return CoSimulationBaseSolver(cosim_solver_settings, level)

class CoSimulationBaseSolver(object):
    """The base class for the CoSimulation Solvers
    The intention is that every solver that derives from this class
    can be used standalone.
    """
    def __init__(self, cosim_solver_settings, level):
        """Constructor of the Base-Solver
        Deriving classes should call it in their constructors
        """
        self.cosim_solver_settings = cosim_solver_settings
        self.lvl = level

        self.io = None
        if "io_settings" in self.cosim_solver_settings:
            self.io = io_factory.CreateIO(self.cosim_solver_settings["io_settings"])

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def AdvanceInTime(self, current_time):
        return current_time + self.cosim_solver_settings["time_step"] # needed if this solver is used as dummy

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

    def ImportData(self, data_name, from_client):
        if self.io:
            self.io.ImportData(data_name, from_client)
    def ImportMesh(self, mesh_name, from_client):
        if self.io:
            self.io.ImportMesh(mesh_name, from_client)

    def ExportData(self, data_name, to_client):
        if self.io:
            self.io.ExportData(data_name, to_client)
    def ExportMesh(self, mesh_name, to_client):
        if self.io:
            self.io.ExportMesh(mesh_name, to_client)

    def MakeDataAvailable(self, data_name, to_client):
        if self.io:
            self.io.MakeDataAvailable(data_name, to_client)
    def MakeMeshAvailable(self, mesh_name, to_client):
        if self.io:
            self.io.MakeMeshAvailable(mesh_name, to_client)

    def GetDataDefinition(self, data_name):
        return self.cosim_solver_settings["data"][data_name]