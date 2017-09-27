from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# setting the domain size for the problem to be solved
domain_size = 2

import math

# importing applications
from KratosMultiphysics import *
from KratosMultiphysics.ErodibleBedApplication import *

# from now on the order is not anymore crucial
#
#


# defining a model part
model_part = ModelPart("SoilPart")

# importing the solver files and adding the variables
import bed_convection_solver
bed_convection_solver.AddVariables(model_part)

# reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("test_out", gid_mode, multifile, deformed_mesh_flag, write_conditions)

model_part_io = ModelPartIO("test")
model_part_io.ReadModelPart(model_part)

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

# importing the solver files
bed_convection_solver.AddDofs(model_part)


for node in model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY_X, 0, 1.0)
    node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
    node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)
    if node.X<2.5:
        node.SetSolutionStepValue(HEIGHT, 0, 1.0)
    else:
        node.SetSolutionStepValue(HEIGHT, 0, 0.0)
    if node.X<0.0001:
        node.Fix(HEIGHT)



# convection solver: CONSTRUCTOR ##########
convection_solver = bed_convection_solver.BedConvectionSolver(model_part, domain_size)
convection_solver.Initialize()

mesh_name=0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(mesh_name, model_part.GetMesh())

gid_io.WriteNodalResults(HEIGHT, model_part.Nodes, 0.0, 0)


delta_t = 0.002
max_time = 3.0

time_to_print=0.0
time=0.0
step=0

while time < max_time:

    time += delta_t
    time_to_print +=delta_t
    print("Current time = ", time)
    model_part.CloneTimeStep(time)

    if(step >= 1):
        convection_solver.Solve()

    # print the results
    mesh_name = time  # if we want the mesh to change at each time step then ****mesh_name = time****

    step = step + 1


    if(time_to_print >= 0.01):
        print("output")
        gid_io.WriteNodalResults(PROJECTED_HEIGHT, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(HEIGHT, model_part.Nodes, time, 0)


        gid_io.Flush()


        time_to_print = 0.0


gid_io.FinalizeResults()
