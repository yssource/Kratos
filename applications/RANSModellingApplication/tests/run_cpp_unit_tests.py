from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # TESTS_OUTPUTS
    # Tester.RunTestSuite("StabilizedAdjointUtilsTestSuite")
    Tester.RunTestSuite("RANSYPlusModels")
    # Tester.RunTestSuite("RANSEvModelsKEpsilonNodalMatrices")
    Tester.RunTestSuite("RANSEvModelsKEpsilonGaussMatrices")
    Tester.RunTestSuite("RANSEvModelsKEpsilonElementResidualMatrices")

if __name__ == '__main__':
    run()