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
        "cad_filename"                  : "tripod_coarse.iga",
        "fem_filename"                  : "tripod.mdpa",
        "fe_refinement_level"           : 0
    },
    "conditions" :
    {
        "general" :
        {
            "apply_integral_method" : false
        },
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
                "exclusive_face_list" : [],
                "penalty_factor"      : 1e1
            }
        },
        "edges" :
        {
            "fe_based" :
            {
                "apply_enforcement_conditions"        : false,
                "penalty_factor_position_enforcement" : 1e1,
                "penalty_factor_tangent_enforcement"  : 1e2,
                "apply_corner_enforcement_conditions" : false,
                "penalty_factor_corner_enforcement"   : 1e4
            },
            "coupling" :
            {
                "apply_coupling_conditions"            : true,
                "penalty_factor_displacement_coupling" : 1e4,
                "penalty_factor_rotation_coupling"     : 1e4
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
        "test_solution" : false
    },
    "regularization" :
    {
        "alpha"             : 1.0,
        "beta"              : 0.001,
        "include_all_poles" : false
    },
    "refinement" :
    {
        "a_posteriori" :
        {
            "apply_a_posteriori_refinement" : true,
            "max_levels_of_refinement"      : 8,
            "mininimum_knot_distance"       : 0.8,
            "fe_point_distance_tolerance"   : 0.1,
            "disp_coupling_tolerance"       : 0.005,
            "rot_coupling_tolerance"        : 0.2
        },
        "a_priori" :
        {
            "apply_a_priori_refinement"         : true,
            "max_levels_of_refinement"          : 10,
            "min_knot_distance_at_max_gradient" : 1.0,
            "exponent"                          : 2
        }
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
print("> Start reconstruction...")
print("========================================================================================================")

start_time = time.time()

cad_mapper = CADMapper(fe_model, cad_model, parameters)
cad_mapper.RunMappingProcess()

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")