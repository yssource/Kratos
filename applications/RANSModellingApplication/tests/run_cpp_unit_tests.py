from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.FAILED_TESTS_OUTPUTS) # TESTS_OUTPUTS
    Tester.RunTestSuite("StabilizedAdjointUtilsTestSuite")
    Tester.RunTestSuite("RANSYPlusModels")
    Tester.RunTestSuite("RANSEvModels")

if __name__ == '__main__':
    run()