from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.EmpireApplication.co_simulation_steady_analysis as CoSimulationSteadyAnalysis

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

import os

# Check other applications dependency
hdf5_is_available = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class TestFSIPotentialSolver(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def TestFSIPotentialSolver(self):
        if not hdf5_is_available:
            self.skipTest("Missing required application: HDF5Application")
        settings_file_name_structure = "ProjectParametersSDoF.json"
        settings_file_name_fluid = "project_parameters_cosim_naca0012_small_fsi.json"
        settings_file_name_coupler = "naca0012_small_parameters_coupling.json"
        work_folder = "fsi_potential_flow_sdof"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name_fluid)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.4919547597104221, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.15880211309589895, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[KratosMultiphysics.REACTION[1]], 30.13222903226336, 0.0, 1e-9)

            for settings_file_name_fluid in os.listdir(os.getcwd()):
                if settings_file_name_fluid.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(settings_file_name_fluid)

    def _runTest(self,settings_file_name):
        model = KratosMultiphysics.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = KratosMultiphysics.Parameters(settings_file.read())


        co_steady_analysis = CoSimulationSteadyAnalysis(model, settings)
        co_steady_analysis.Run()
        self.main_model_part = model.GetModelPart(settings["solver_settings"]["model_part_name"].GetString())

    def _check_results(self, result, reference, rel_tol, abs_tol):
        isclosethis = t_isclose(result, reference, rel_tol, abs_tol)

        full_msg =  "Failed with following parameters:\n"
        full_msg += str(result) + " != " + str(reference) + ", rel_tol = "
        full_msg += str(rel_tol) + ", abs_tol = " + str(abs_tol)

        self.assertTrue(isclosethis, msg=full_msg)

if __name__ == '__main__':
    UnitTest.main()