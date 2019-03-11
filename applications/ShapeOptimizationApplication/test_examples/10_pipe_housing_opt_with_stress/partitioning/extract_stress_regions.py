import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosCSM
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

import time, math

#############################################################################
# Read mdpa
#############################################################################

model = Kratos.Model()

csm_part_name = "CSM_domain_with_large_hole_and_cutout_and_bolts"
csm_part_filename = csm_part_name
csm_part = model.CreateModelPart(csm_part_name)
# csm_part.AddNodalSolutionStepVariable(KratosShape.VECTOR_VARIABLE)
# csm_part.AddNodalSolutionStepVariable(KratosShape.VECTOR_VARIABLE_MAPPED)
# csm_part.AddNodalSolutionStepVariable(KratosShape.SHAPE_UPDATE)

model_part_io = Kratos.ModelPartIO(csm_part_filename)
model_part_io.ReadModelPart(csm_part)
csm_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE,3)
csm_part.Properties[0].SetValue(Kratos.CONSTITUTIVE_LAW, KratosCSM.LinearElastic3DLaw())
csm_part.SetBufferSize(1)

#############################################################################
# Identify and partition different stress areas
#############################################################################

print("> Starting to partition stress regions...")
start_time = time.time()

# top_lower_elements = []
# top_upper_elements = []
# bottom_lower_elements = []
# bottom_top_elements = []

top_lower_index = 1
top_upper_index = 2
bottom_lower_index = 3
bottom_upper_index = 4

others_index = 0

# # Clear partition index
# for elem in csm_part.Elements:
#     elem.SetValue(Kratos.PARTITION_INDEX,0)

outer_radius = 0.08
inner_radius = 0.05

# First identify all top elements
for elem in csm_part.Elements:

    elem_center = [0,0,0]
    for node in elem.GetNodes():
        elem_center[0] += node.X
        elem_center[1] += node.Y
        elem_center[2] += node.Z
    num_nodes = len(elem.GetNodes())
    elem_center[0] /= num_nodes
    elem_center[1] /= num_nodes
    elem_center[2] /= num_nodes

    is_element_in_top_region = False
    if elem_center[0] > 0:
        is_element_in_top_region = True

    if is_element_in_top_region:

        if elem_center[2] > 0.11 and elem_center[2] < 0.16:
            if elem_center[0] < outer_radius:
                if elem_center[0] < inner_radius:
                    elem.SetValue(Kratos.PARTITION_INDEX,top_lower_index)
                else:
                    elem.SetValue(Kratos.PARTITION_INDEX,top_upper_index)

        elif elem_center[2] > 0.16:
            r_actual_square = (elem_center[0]-0)**2 + (elem_center[2]-0.16)**2

            if r_actual_square < outer_radius**2:
                if r_actual_square < inner_radius**2:
                    elem.SetValue(Kratos.PARTITION_INDEX,top_lower_index)
                else:
                    elem.SetValue(Kratos.PARTITION_INDEX,top_upper_index)
    else:

        if elem_center[2] > 0.11 and elem_center[2] < 0.16:
            if elem_center[0] > -outer_radius:
                if elem_center[0] > -inner_radius:
                    elem.SetValue(Kratos.PARTITION_INDEX,bottom_lower_index)
                else:
                    elem.SetValue(Kratos.PARTITION_INDEX,bottom_upper_index)

        elif elem_center[2] > 0.16:
            r_actual_square = (elem_center[0]-0)**2 + (elem_center[2]-0.16)**2

            if r_actual_square < outer_radius**2:
                if r_actual_square < inner_radius**2:
                    elem.SetValue(Kratos.PARTITION_INDEX,bottom_lower_index)
                else:
                    elem.SetValue(Kratos.PARTITION_INDEX,bottom_upper_index)

print("> Finished partitioning in" ,round( time.time()-start_time, 3 ), " s.")

#############################################################################
# Output results
#############################################################################

print("> Starting to output partitions...")
start_time = time.time()

# Write SubModelParts for each partition index
with open('partitions.txt','w') as file:
    for itr in range(4):
        file.write("Begin SubModelPart stress_partition_"+str(itr+1)+"\n")
        file.write("    Begin SubModelPartElements\n")
        for elem in csm_part.Elements:
            if elem.GetValue(Kratos.PARTITION_INDEX) == itr+1:
                elem_id = str(elem.Id)
                file.write("        "+elem_id+"\n")
        file.write("    End SubModelPartElements\n")
        file.write("End SubModelPart\n\n")

vtk_parameters = Kratos.Parameters("""
{
    "model_part_name": "Parts_structure",
    "file_format": "ascii",
    "output_precision": 7,
    "output_control_type": "step",
    "output_frequency": 1.0,
    "output_sub_model_parts": true,
    "folder_name": "partition_results",
    "save_output_files_in_folder": true,
    "nodal_solution_step_data_variables": [],
    "nodal_data_value_variables": [],
    "element_data_value_variables": ["PARTITION_INDEX"],
    "condition_data_value_variables": []
}""")

vtk_output = VtkOutputProcess(model, vtk_parameters)
vtk_output.ExecuteInitialize()
vtk_output.ExecuteBeforeSolutionLoop()
vtk_output.ExecuteInitializeSolutionStep()
vtk_output.PrintOutput()
vtk_output.ExecuteFinalizeSolutionStep()
vtk_output.ExecuteFinalize()

print("> Finished output of partitions in" ,round( time.time()-start_time, 3 ), " s.")