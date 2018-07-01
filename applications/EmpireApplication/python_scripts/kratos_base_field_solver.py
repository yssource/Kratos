from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
import io_factory

class KratosBaseFieldSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings):
        super(KratosBaseFieldSolver, self).__init__(cosim_solver_settings)

        self.model = KratosMultiphysics.Model()

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.project_parameters = KratosMultiphysics.Parameters(parameters)

    def Initialize(self):
        self._GetAnalysisStage().Initialize()

    def Finalize(self):
        self._GetAnalysisStage().Finalize()

    def AdvanceInTime(self, current_time):
        return self._GetAnalysisStage()._GetSolver().AdvanceInTime(current_time)

    def Predict(self):
        self._GetAnalysisStage()._GetSolver().Predict()

    def InitializeSolutionStep(self):
        self._GetAnalysisStage().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetAnalysisStage().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetAnalysisStage().OutputSolutionStep()

    def SolveSolutionStep(self):
        self._GetAnalysisStage()._GetSolver().SolveSolutionStep()

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

    def _GetAnalysisStage(self):
        if not hasattr(self, '_analysis_stage'):
            self._analysis_stage = self._CreateAnalysisStage()
        return self._analysis_stage

    def _CreateAnalysisStage(self):
        raise Exception("Creation of the AnalysisStage must be implemented in the derived class!")