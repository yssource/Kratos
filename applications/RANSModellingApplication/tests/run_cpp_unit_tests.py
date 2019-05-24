from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.FAILED_TESTS_OUTPUTS) # TESTS_OUTPUTS
    # Tester.RunTestSuite("RansCalculationUtilities")
    # Tester.RunTestSuite("RANSYPlusModels")
    # Tester.RunTestSuite("RANSEvModelsKEpsilonNodalMatrices")
    # Tester.RunTestSuite("RANSEvModelsKEpsilonGaussMatrices")
    Tester.RunTestSuite("RANSEvModelsKEpsilonElementResidualMatrices")
    Tester.RunTestSuite("RANSEvModelsKEpsilonElementResidualMatricesT1")

if __name__ == '__main__':
    run()