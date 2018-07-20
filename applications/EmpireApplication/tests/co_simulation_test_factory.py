from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EmpireApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import os

import co_simulation_test_case

class TestSmallCoSimulationCases(KratosUnittest.TestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''
    def test_MokFSI(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

class TestCoSimulationCases(KratosUnittest.TestCase):
    '''This class contains "full" CoSimulation-Cases, too large for the nightly suite and therefore
    have to be in the validation-suite
    '''
    def test_WallFSI(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

if __name__ == '__main__':
    KratosUnittest.main()
