import sys

print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))

from KratosMultiphysics import *
Tester.DisableAllTestCases() # disabling the tests of the core

application_module_name = "KratosMultiphysics." + sys.argv[1]
__import__(application_module_name)
Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
print(Tester.RunAllTestCases())
