# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosCSM
import math

################################################################################################
# Read fe data
################################################################################################
fe_model = KratosMultiphysics.Model()

fe_model_part = fe_model.CreateModelPart("origin_part")
fe_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LAMBDA)
fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

fem_input_filename = "3D_shell_in_3D_space"
model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename)
model_part_io.ReadModelPart(fe_model_part)

################################################################################################
# Create some scalar field lambda
################################################################################################
for node in fe_model_part.Nodes:

    amplitude = 25
    inner_parabula_radius = 5
    outer_parabula_radius = 40

    a = amplitude/((inner_parabula_radius-outer_parabula_radius)**2)
    b = -2*a*outer_parabula_radius
    c = a*outer_parabula_radius**2

    x = math.sqrt(node.X**2 + node.Y**2 + node.Z**2)

    if x <= outer_parabula_radius:
        height = a*x**2 + b*x + c
    else:
        height = 0

    node.SetSolutionStepValue(KratosMultiphysics.LAMBDA, height)

################################################################################################
# Compute gradient of displacement
################################################################################################
local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(fe_model_part, KratosMultiphysics.LAMBDA, KratosMultiphysics.SHAPE_SENSITIVITY, KratosMultiphysics.NODAL_AREA)
local_gradient.Execute()

################################################################################################
# Output results
################################################################################################
nodal_variables = ["LAMBDA", "SHAPE_SENSITIVITY", "NODAL_AREA"]

from gid_output import GiDOutput
nodal_results=nodal_variables
gauss_points_results=[]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = False
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

gig_io = GiDOutput("results", VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
gig_io.initialize_results(fe_model_part)
gig_io.write_results(1, fe_model_part, nodal_results, gauss_points_results)
gig_io.finalize_results()