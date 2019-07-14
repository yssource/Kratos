# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.MGISApplication as MGISApplication
#import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests or test_classes to create the suits

##### SELF-CONTAINED TESTS #####
# Processes test

##### SMALL TESTS #####

##### NIGHTLY TESTS #####

##### VALIDATION TESTS #####

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    # These tests are executed by the continuous integration tool, so they have to be very fast!
    # Execution time << 1 sec on a regular PC !!!
    # If the tests in the smallSuite take too long then merging to master will not be possible!
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    nightSuite = suites['nightly'] # These tests are executed in the nightly build

    ### Adding the self-contained tests
    # Constitutive Law tests
    #smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestConstitutiveLaw]))

    ### Adding Small Tests
    nightSuite.addTests(smallSuite)

    ### Adding Validation Tests
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    validationSuite.addTests(allSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
