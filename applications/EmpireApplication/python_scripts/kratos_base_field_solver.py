from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from co_simulation_base_solver import CoSimulationBaseSolver


class KratosBaseFieldSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(KratosBaseFieldSolver, self).__init__(cosim_solver_settings, level)
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
        new_time = self._GetAnalysisStage()._GetSolver().AdvanceInTime(current_time)
        self._GetAnalysisStage().time = new_time # only needed to print the time correctly
        return new_time

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


    def _GetAnalysisStage(self):
        if not hasattr(self, '_analysis_stage'):
            self._analysis_stage = self._CreateAnalysisStage()
        return self._analysis_stage

    def _CreateAnalysisStage(self):
        raise Exception("Creation of the AnalysisStage must be implemented in the derived class!")