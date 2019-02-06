# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape

# Additional imports
import time
import os
from cad_reconstruction import CADMapper

# Parameters
parameters = KratosMultiphysics.Parameters("""
{
    "input" :
    {
        "cad_filename"                  : "plate_fine_embedded.iga",
        "fem_filename"                  : "plate.mdpa",
        "fe_refinement_level"           : 0,
        "variable_to_map"               : "SHAPE_CHANGE"
    },
    "conditions" :
    {
        "apply_integral_method" : false,
        "faces" :
        {
            "curvature" :
            {
                "apply_curvature_minimization" : false,
                "penalty_factor"               : 1e-1
            },
            "mechanical" :
            {
                "apply_KL_shell"      : false,
                "exclusive_face_list" : ["Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<7>"],
                "penalty_factor"      : 1e3
            },
            "rigid" :
            {
                "apply_rigid_conditions" : false,
                "exclusive_face_list"    : ["Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<5>.Hole"],
                "penalty_factor"         : 1e5
            }
        },
        "edges" :
        {
            "fe_based" :
            {
                "apply_enforcement_conditions"        : true,
                "penalty_factor_tangent_enforcement"  : 1000,
                "penalty_factor_position_enforcement" : 1000,
                "apply_corner_enforcement_conditions" : true,
                "penalty_factor_corner_enforcement"   : 10000
            },
            "coupling" :
            {
                "apply_coupling_conditions" : false
            }
        }
    },
    "point_search" :
    {
        "boundary_tessellation_tolerance" : 0.01,
        "patch_bounding_box_tolerance"    : 1.0
    },
    "solution" :
    {
        "iterations"    : 1,
        "test_solution" : false
    },
    "regularization" :
    {
        "alpha" : 0.1,
        "beta"  : 0.001
    },
    "output":
    {
        "results_directory"           : "01_Results",
        "resulting_geometry_filename" : "reconstructed_geometry.iga"
    }
}""")

fe_model = KratosMultiphysics.Model()
cad_model = None

print("\n\n========================================================================================================")
print("> Start some preprocessing...")
print("========================================================================================================")

# Read FE data
fem_input_filename = parameters["input"]["fem_filename"].GetString()
name_variable_to_map = parameters["input"]["variable_to_map"].GetString()
variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

fe_model_part = fe_model.CreateModelPart("origin_part")
fe_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
fe_model_part.AddNodalSolutionStepVariable(variable_to_map)
fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)

model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
model_part_io.ReadModelPart(fe_model_part)

# Create some displacement field
for node in fe_model_part.GetSubModelPart("all_faces").Nodes:
    import math

    # amplitude = 25
    # inner_ring_radius = 20
    # sigma = inner_ring_radius/2
    # gauss_value = amplitude*math.exp(-(node.X**2/(2*sigma**2)+node.Y**2/(2*sigma**2)))
    # node.SetSolutionStepValue(variable_to_map, [0,0,gauss_value])

    amplitude = 25
    inner_parabula_radius = 5
    outer_parabula_radius = 40
    # outer_parabula_radius = 30 # Radius to center of holes

    a = amplitude/((inner_parabula_radius-outer_parabula_radius)**2) # double a to have a wider parabula
    b = -2*a*outer_parabula_radius
    c = a*outer_parabula_radius**2

    x = math.sqrt(node.X**2 + node.Y**2)

    if x <= outer_parabula_radius:
        y = a*x**2 + b*x + c
    else:
        y = 0

    node.SetSolutionStepValue(variable_to_map, [0,0,y])


print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

start_time = time.time()

cad_mapper = CADMapper(fe_model, cad_model, parameters)

cad_mapper.Initialize()
cad_mapper.Map()
cad_mapper.Finalize()

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")