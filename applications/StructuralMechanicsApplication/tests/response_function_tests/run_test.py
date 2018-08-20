# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_response_function_factory

import KratosMultiphysics.kratos_utilities as kratos_utils


def _get_test_working_dir():
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    return this_file_dir

class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

def assertAlmostEqual(a, b, tol=7):
    if abs(a-b) > 10**(-tol):
        raise RuntimeError("ERROR: "+str(a)+"!="+str(b))

def assertEqualSignificantDigits(a, b, digits=7):
    epsilon = 10 ** (-digits)
    if b == 0.0:
        if abs(a) > epsilon:
            raise RuntimeError("ERROR: "+str(a)+"!="+str(b))
    elif a == 0.0:
        if abs(b) > epsilon:
            raise RuntimeError("ERROR: "+str(a)+"!="+str(b))
    elif abs(a/b - 1) > epsilon:
        raise RuntimeError("ERROR: "+str(a)+"!="+str(b))

def check_primal_results(model_part):
    assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X), 0.0)
    assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
    assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0023682621559198303)
    assertAlmostEqual(model_part.Nodes[4].GetSolutionStepValue(KratosMultiphysics.ROTATION_Y), -0.005700479161950828)

class StructuralResponseFunctionTestFactory():

    file_name = None

    def test_execution(self):
        with controlledExecutionScope(_get_test_working_dir()):
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters( parameter_file.read())

            # To avoid many prints
            if (parameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            problem_name = parameters["problem_data"]["problem_name"].GetString()
            model_part = KratosMultiphysics.ModelPart(problem_name)

            self.response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model_part)

            # import model part
            model_part_io = KratosMultiphysics.ModelPartIO(parameters["problem_data"]["problem_name"].GetString())
            model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
            model_part_io.ReadModelPart(model_part)

            # call response function
            self.response_function.Initialize()

            self.response_function.InitializeSolutionStep()

            # check results of primal analysis on primal model part
            check_primal_results(self.response_function.primal_analysis.model.GetModelPart("rectangular_plate_structure"))
            # check primal results on adjoint model part
            check_primal_results(self.response_function.adjoint_analysis.model.GetModelPart("rectangular_plate_structure"))

            self.response_function.CalculateValue()
            self.value = self.response_function.GetValue()
            self.response_function.CalculateGradient()

            # get the final results of the sensitivity analysis
            self.gradient = self.response_function.GetShapeGradient()

            self.check_adjoint_results()

            self.response_function.FinalizeSolutionStep()

            self.response_function.Finalize()

    def tearDown(self):
        with controlledExecutionScope(_get_test_working_dir()):
            kratos_utils.DeleteFileIfExisting("rectangular_plate_structure.post.bin")
            kratos_utils.DeleteFileIfExisting("rectangular_plate_structure.time")
            kratos_utils.DeleteFileIfExisting("rectangular_plate_structure.h5")
            kratos_utils.DeleteFileIfExisting("rectangular_plate_structure-1.0000.h5")
            kratos_utils.DeleteFileIfExisting("response_function_tests.post.lst")

class TestAdjointStrainEnergyResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_strain_energy_response"

    def check_adjoint_results(self):
        assertAlmostEqual(self.value, 0.011841310779599155)

        nodeId = 4
        assertAlmostEqual(self.gradient[nodeId][0], -0.02211883892695645)
        assertAlmostEqual(self.gradient[nodeId][1], 0.0034659667625643337)
        assertAlmostEqual(self.gradient[nodeId][2], -3.942528133752394e-07, 12)

class TestAdjointDisplacementResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_displacement_response"
    def check_adjoint_results(self):
        model_part = self.response_function.adjoint_analysis.model.GetModelPart("rectangular_plate_structure")
        assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_X), 0.00011007358655353125)
        assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y), -0.00019567412134447555)
        assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Z), -3216405525244.02, 0)
        assertAlmostEqual(model_part.Nodes[4].GetSolutionStepValue(ADJOINT_ROTATION_Y), 8198364044238.647, 1)

        assertAlmostEqual(self.value, 0.0023682621559198303)

        nodeId = 4
        assertAlmostEqual(self.gradient[nodeId][0], 35162269878424.273, 1)
        assertAlmostEqual(self.gradient[nodeId][1], 4507732710819.046, 1)
        assertAlmostEqual(self.gradient[nodeId][2], 1125096436.1808963, 3)

class TestAdjointStressResponseFunction(StructuralResponseFunctionTestFactory):
    file_name = "adjoint_stress_response"
    def check_adjoint_results(self):
        model_part = self.response_function.adjoint_analysis.model.GetModelPart("rectangular_plate_structure")
        assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_X), -3.4423476111614904e-17, 12)
        assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y), 6.278340199159814e-17, 12)
        assertAlmostEqual(model_part.Nodes[5].GetSolutionStepValue(ADJOINT_DISPLACEMENT_Z), 0.9126590318076583)
        assertAlmostEqual(model_part.Nodes[4].GetSolutionStepValue(ADJOINT_ROTATION_Y), -1.6105506659693518)

    def
        assertAlmostEqual(self.value, -0.8233392989483465)

        nodeId = 4
        assertAlmostEqual(self.gradient[nodeId][0], -8.480870212810133)
        assertAlmostEqual(self.gradient[nodeId][1], -2.6870951672193213)
        assertAlmostEqual(self.gradient[nodeId][2], -0.0002913758359272386, 9)


if __name__ == "__main__":
    results = []
    tests = [TestAdjointStressResponseFunction(),TestAdjointDisplacementResponseFunction(),TestAdjointStrainEnergyResponseFunction()]
    for t in tests:
        try:
            t.test_execution()
        except Exception as e:
            results.append([t.__class__.__name__, e])
        else:
            results.append([t.__class__.__name__, "SUCCESS"])
        finally:
            t.tearDown()

    for r in results:
        print(r)