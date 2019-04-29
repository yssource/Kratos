from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_LIST) # TESTS_OUTPUTS
    Tester.RunTestSuite("StabilizedAdjointUtilsTestSuite")
    Tester.RunTestSuite("RANSYPlusModels")

if __name__ == '__main__':
    run()