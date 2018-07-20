from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.EmpireApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import os

import co_simulation_test_case

class TestKratosSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_KratosStructuralMechanicsSolver(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")

    def test_KratosFluidDynamicsSolver(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")


class TestSDofSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_SDofSolver(self):
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")


class TestEmpireSolver(co_simulation_test_case.CoSimulationTestCase):
    def test_EmpireSolverWrapper(self):
        if "EMPIRE_API_LIBSO_ON_MACHINE" not in os.environ:
            self.skipTest("EMPIRE is not available")
        with co_simulation_test_case.ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # self.createTest('test_structural_mesh_motion_2d/rectangle_2D3N_test')
            # self.runTest()
            kratos_utils.DeleteFileIfExisting("./test_mdpa_files/rectangle_2D3N_test.time")


if __name__ == '__main__':
    KratosUnittest.main()
