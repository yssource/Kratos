# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from gid_output_process import GiDOutputProcess
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import math

#############################################################################
# Parameters
#############################################################################

box_mdpa_filename = "box"
sphere_mdpa_filename = "sphere"
skin_mdpa_filename = "skin"
levels_of_refinement = 3
sphere_radius = 0.4

output_parameters = Parameters("""
{
    "result_file_configuration" : {
        "gidpost_flags"       : {
            "GiDPostMode"           : "GiD_PostBinary",
            "WriteDeformedMeshFlag" : "WriteDeformed",
            "WriteConditionsFlag"   : "WriteConditions",
            "MultiFileFlag"         : "SingleFile"
        },
        "file_label"          : "step",
        "output_control_type" : "step",
        "output_frequency"    : 1,
        "body_output"         : true,
        "node_output"         : false,
        "skin_output"         : false,
        "nodal_results"       : [],
        "gauss_point_results" : []
    },
    "point_data_configuration"  : []
}""")

#############################################################################
# Define general functions
#############################################################################
def ReadMdpa( mdpa_filename ):
    mdpa = ModelPart(mdpa_filename)
    mdpa.AddNodalSolutionStepVariable(DISTANCE)
    mdpa.AddNodalSolutionStepVariable(VELOCITY)
    model_part_io = ModelPartIO(mdpa_filename)
    model_part_io.ReadModelPart(mdpa)
    return mdpa

def OutputMdpa( output_mdpa, output_filename, parameters ):
    gid_output_original = GiDOutputProcess( output_mdpa, output_filename, parameters )
    gid_output_original.ExecuteInitialize()
    gid_output_original.ExecuteBeforeSolutionLoop()
    gid_output_original.ExecuteInitializeSolutionStep()
    gid_output_original.PrintOutput()
    gid_output_original.ExecuteFinalizeSolutionStep()
    gid_output_original.ExecuteFinalize()

def EvaluateImplicitRepresentation(skin_mdpa):
    error_vector = []
    for node in skin_mdpa.Nodes:
        error = math.sqrt(node.X**2 + node.Y**2 + node.Z**2) - sphere_radius
        error_vector.append(error)

    norm_2 = 0.0
    for entry in error_vector:
        norm_2 = norm_2 + entry*entry
    norm_2 = math.sqrt(norm_2)

    return norm_2

#############################################################################
# Main body
#############################################################################

number_of_elements = []
error_values = []

# Defining the model_parts
box_mdpa = ReadMdpa(box_mdpa_filename)
sphere_mdpa = ReadMdpa(sphere_mdpa_filename)

box_mdpa.Properties[0].SetValue(CONSTITUTIVE_LAW, LinearElastic3DLaw())
box_mdpa.SetBufferSize(2)

# Refine
for refinement_level in range(0,levels_of_refinement):

    print("#########################")
    print("> Starting refinement level ", refinement_level)

    # Calculate signed distance process
    print("> Calculating distances!!!!!")
    calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(sphere_mdpa, box_mdpa)
    calculate_distance_process.Execute()
    # OutputMdpa(box_mdpa, "box_after_distance_calculation", output_parameters)

    # Generate reconstruction
    print("> Generating skin mdpa!!!!!")
    skin_mdpa = ModelPart(skin_mdpa_filename)
    calculate_distance_process.GenerateSkinModelPart(skin_mdpa)
    # OutputMdpa(skin_mdpa, skin_mdpa_filename, output_parameters)

    # Evaluate error
    number_of_elements.append(box_mdpa.NumberOfElements())
    error_values.append(EvaluateImplicitRepresentation(skin_mdpa))

    # Mesh refinement requires identification of elemental and nodal neighbours
    FindNodalNeighboursProcess(box_mdpa, 100, 100).Execute()
    FindElementalNeighboursProcess(box_mdpa,3,100).Execute()

    # For only adaptive refinement: Note that all split elements are already marked as split using the variable SPLIT_ELEMENT

    # Refinement of complete box mesh
    for element in box_mdpa.Elements:
        element.SetValue(SPLIT_ELEMENT,True)

    # Actual refinement
    refine_on_reference = False
    interpolate_internal_variables = False

    refinement_tool = LocalRefineTetrahedraMesh(box_mdpa)
    refinement_tool.LocalRefineMesh(refine_on_reference, interpolate_internal_variables)