from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EmpireApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from co_simulation_analysis import CoSimulationAnalysis

import os, json

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class CoSimulationTestCase(KratosUnittest.TestCase):

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            self.cosim_parameters = json.load(parameter_file)

        # # To avoid many prints
        # echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        # if (echo_level == 0):
        #     KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def runTest(self):
        CoSimulationAnalysis(self.cosim_parameters).Run()
