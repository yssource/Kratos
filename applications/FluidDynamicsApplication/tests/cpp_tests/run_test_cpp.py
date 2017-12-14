from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
# from KratosMultiphysics.mpi import *
# from KratosMultiphysics.MPISearchApplication import *

Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)

# Tester.RunTestCases('*GenerateBins3DMPI')

# Tester.RunTestCases('TestTriangle2D3*')
# Tester.RunTestCases('*Embedded*')
# Tester.RunTestCases('*Navier*')
Tester.RunTestCases('*Compressible*')
#Tester.RunTestCases('*CompressibleNavierStokes2D3NConstant*')
# Tester.RunTestCases('*Intersected*')
# Tester.RunTestCases('*Distance*')
# Tester.RunTestCases('*DiscontUtils')
#Tester.RunTestCases('*DivideGeometry*')
#Tester.RunTestCases('*ModifiedShapeFunctions*')

# Tester.RunTestSuite('FluidDynamicsApplicationFastSuite')
# Tester.RunTestSuite('KratosCoreGeometriesFastSuite')
# Tester.RunTestSuite('KratosCoreFastSuite')

# Tester.RunAllTestCases()
