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
    "inpute" :
    {
        "cad_filename"                  : "tripod.iga",
        "fem_filename"                  : "tripod.mdpa",
        "fe_refinement_level"           : 0,
        "variable_to_map"               : "SHAPE_CHANGE"
    },
    "points_projection" :
    {
        "boundary_tessellation_tolerance" : 0.01,
        "patch_bounding_box_tolerance"    : 1.0
    },
    "boundary_conditions" :
    {
        "penalty_factor" : 10
    },
    "solution" :
    {
        "iterations"    : 1,
        "test_solution" : true
    },
    "regularization" :
    {
        "beta" : 0.01
    },
    "output":
    {
        "results_directory"           : "01_Results",
        "resulting_geometry_filename" : "reconstructed_geometry.iga"
    }
}""")


print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

start_time = time.time()

fe_model = KratosMultiphysics.Model()
cad_model = None

cad_mapper = CADMapper(fe_model, cad_model, parameters)

cad_mapper.Initialize()
cad_mapper.Map()
cad_mapper.Finalize()

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")