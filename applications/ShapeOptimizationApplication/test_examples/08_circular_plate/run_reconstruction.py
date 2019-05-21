# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape

# Additional imports
import time
import os
from cad_reconstruction_mapper import CADMapper
import ANurbs as an

# Parameters
parameters = KratosMultiphysics.Parameters("""
{
    "input" :
    {
        "cad_filename"                  : "plate_with_fine_strips.iga",
        "fem_filename"                  : "plate.mdpa",
        "fe_refinement_level"           : 0
    },
    "conditions" :
    {
        "general" :
        {
            "apply_integral_method" : false,
            "mapping_cad_fem"       :
            {
                "mapper_type"   : "nearest_element",
                "search_radius" : 10.0,
                "echo_level"    : 0
            }
        },
        "faces" :
        {
            "mechanical" :
            {
                "apply_KL_shell"      : true,
                "exclusive_face_list" : [
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<11>",
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<9>",
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<4>",
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<6>",
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<7>"
                ],
                "penalty_factor"      : 1e1
            },
            "rigid" :
            {
                "apply_rigid_conditions" : true,
                "exclusive_face_list"    : [
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<5>.Hole",
                    "Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepFace<11>.Hole"
                ],
                "penalty_factor"         : 1e6
            }
        },
        "edges" :
        {
            "direct" :
            {
                "apply_enforcement_conditions" : true,
                "exclusive_edge_list": [
                    ["Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepEdge<0>", [0,0,25]],
                    ["Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepEdge<5>" , [0,0,25]],
                    ["Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepEdge<9>" , [0,0,25]],
                    ["Rhino<5b74e2e3-fd33-446e-a112-0480ea75a27d>.BrepEdge<13>" , [0,0,25]]
                ],
                "penalty_factor_position_enforcement" : 1e3
            },
            "fe_based" :
            {
                "apply_enforcement_conditions"        : false,
                "penalty_factor_position_enforcement" : 1e3,
                "penalty_factor_tangent_enforcement"  : 1e3,
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
        "alpha" : 0.0,
        "beta"  : 0.001
    },
    "refinement" :
    {
        "a_posteriori" :
        {
            "apply_a_posteriori_refinement" : false,
            "max_levels_of_refinement"      : 8,
            "mininimum_knot_distance"       : 0.1,
            "fe_point_distance_tolerance"   : 1.0,
            "disp_coupling_tolerance"       : 0.01,
            "rot_coupling_tolerance"        : 0.4
        },
        "a_priori" :
        {
            "apply_a_priori_refinement"         : false,
            "max_levels_of_refinement"          : 10,
            "min_knot_distance_at_max_gradient" : 10.0,
            "exponent"                          : 2
        }
    },
    "output":
    {
        "results_directory"           : "01_Results",
        "resulting_geometry_filename" : "reconstructed_geometry",
        "echo_level"                  : 5
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
fe_model_part.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE_ABSOLUTE)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.GRAD_SHAPE_CHANGE)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.FITTING_ERROR)
fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)

fem_input_filename = parameters["input"]["fem_filename"].GetString()
model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
model_part_io.ReadModelPart(fe_model_part)

KratosShape.GeometryUtilities(fe_model_part).ComputeUnitSurfaceNormals(True)

# Create some displacement field
for node in fe_model_part.GetSubModelPart("all_faces").Nodes:
    import math

    # parabula on the inside
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
    node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE, [0,0,normal[2]*height])

    # step on the outside
    plate_radius = 90
    transition_start_radius = 60
    step_start_radius = 80
    height = 0.5*amplitude
    if x > transition_start_radius:
        if x > step_start_radius:
            node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE, [0,0,height])
        else:
            local_coord = (x-transition_start_radius)/(step_start_radius-transition_start_radius) * math.pi # varies 0 and pi
            h = (1 - math.cos(local_coord))*height*0.5
            node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE, [0,0,h])

    # plate_radius = 90
    # step_start_radius = 60
    # height = 0.5*amplitude
    # if x > step_start_radius:
    #     local_coord = (x-step_start_radius)/(plate_radius-step_start_radius) * math.pi # varies 0 and pi/4
    #     h = (1 - math.cos(local_coord))*height
    #     node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE, [0,0,h])


print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

start_time = time.time()

cad_mapper = CADMapper(fe_model, cad_model, parameters)
cad_mapper.RunMappingProcess()

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")