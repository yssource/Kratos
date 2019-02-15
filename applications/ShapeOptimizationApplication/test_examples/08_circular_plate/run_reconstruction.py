# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape

# Additional imports
import time
import os
from cad_reconstruction import CADMapper
import ANurbs as an

# Parameters
parameters = KratosMultiphysics.Parameters("""
{
    "input" :
    {
        "cad_filename"                  : "plate_fine_rotated.iga",
        "fem_filename"                  : "plate_rotated.mdpa",
        "fe_refinement_level"           : 0
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
                "exclusive_face_list" : ["Rhino<9f4d24e9-d364-41c2-8a1b-a28ccee54f7b>.BrepFace<0>"],
                "penalty_factor"      : 1e3
            },
            "rigid" :
            {
                "apply_rigid_conditions" : false,
                "exclusive_face_list"    : ["Rhino<9f4d24e9-d364-41c2-8a1b-a28ccee54f7b>.BrepFace<2>.Hole"],
                "penalty_factor"         : 1e5
            }
        },
        "edges" :
        {
            "fe_based" :
            {
                "apply_enforcement_conditions"        : false,
                "penalty_factor_position_enforcement" : 1e3,
                "penalty_factor_tangent_enforcement"  : 0.0,
                "apply_corner_enforcement_conditions" : false,
                "penalty_factor_corner_enforcement"   : 1e4
            },
            "coupling" :
            {
                "apply_coupling_conditions"            : true,
                "penalty_factor_displacement_coupling" : 1e3,
                "penalty_factor_rotation_coupling"     : 1e3
            }
        }
    },
    "drawing_parameters" :
    {
        "cad_drawing_tolerance"           : 1e-3,
        "boundary_tessellation_tolerance" : 1e-2,
        "patch_bounding_box_tolerance"    : 1.0,
        "min_span_length"                 : 1e-7
    },
    "solution" :
    {
        "iterations"    : 1,
        "test_solution" : true
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
cad_model = an.Model()

print("\n\n========================================================================================================")
print("> Start some preprocessing...")
print("========================================================================================================")

# Read FE data
fe_model_part = fe_model.CreateModelPart("origin_part")
fe_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)
fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)

fem_input_filename = parameters["input"]["fem_filename"].GetString()
model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
model_part_io.ReadModelPart(fe_model_part)

KratosShape.GeometryUtilities(fe_model_part).ComputeUnitSurfaceNormals(True)

# Create some displacement field
for node in fe_model_part.GetSubModelPart("all_faces").Nodes:
    import math

    # amplitude = 25
    # inner_ring_radius = 20
    # sigma = inner_ring_radius/2
    # gauss_value = amplitude*math.exp(-(node.X**2/(2*sigma**2)+node.Y**2/(2*sigma**2)))
    # node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE, [0,0,gauss_value])

    amplitude = 25
    inner_parabula_radius = 5
    outer_parabula_radius = 40
    # outer_parabula_radius = 30 # Radius to center of holes

    a = amplitude/((inner_parabula_radius-outer_parabula_radius)**2) # double a to have a wider parabula
    b = -2*a*outer_parabula_radius
    c = a*outer_parabula_radius**2

    x = math.sqrt(node.X**2 + node.Y**2 + node.Z**2)

    if x <= outer_parabula_radius:
        height = a*x**2 + b*x + c
    else:
        height = 0

    normal = node.GetSolutionStepValue(KratosShape.NORMALIZED_SURFACE_NORMAL)

    node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE, [normal[0]*height,normal[1]*height,normal[2]*height])


print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

start_time = time.time()

cad_mapper = CADMapper(fe_model, cad_model, parameters)
cad_mapper.RunMappingProcess()

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")