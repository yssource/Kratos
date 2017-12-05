# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import time

# ======================================================================================================================================
# Parameters
# ====================================================================================================================================== 

# Parameters
reconstruction_parameters = Parameters("""
{
    "inpute_parameters": 
    {
        "cad_geometry_filename"         : "tripod_geometry.json",
        "cad_integration_data_filename" : "tripod_integration_data.json",
        "fem_filename"                  : "tripod",
        "fe_refinement_level"           : 0,
        "shape_change_variable_name"    : "SHAPE_CHANGE_ABSOLUTE"
    },
    "solution_parameters" : 
    {
        "strategy"            : "displacement_mapping",
        "strategy_specifc_parameters": 
        {
            "fem_gauss_integration_degree" : 5
        },
        "solution_iterations" : 1,
        "constraints"         : 
        {
            "set_displacement_coupling_on_all_coupling_points"  : true,
            "penalty_factor_for_displacement_coupling"          : 1e3,
            "set_rotation_coupling_on_all_coupling_points"      : true,            
            "penalty_factor_for_rotation_coupling"              : 1e3,
            "set_dirichlet_constraints"                         : false,
            "list_of_edge_ids_with_dirichlet_constraints"       : [],                               
            "penalty_factor_for_dirichlet_constraints"          : 1e3,
            "set_constraint_to_enforce_tangent_continuity"      : false,
            "list_of_edge_ids_with_tangent_constraints"         : [],                                           
            "penalty_factor_for_tangent_continuity_constraints" : 1e3,
            "penalty_multiplier"                                : 1.0         
        },
        "projection_parameters": 
        {
            "projection_strategy"                               : "single_search_tree",
            "search_radius_with_multiple_trees"                 : 1.0,               
            "automatic_initialization_using_greville_abscissae" : true,
            "refinement_iterations_of_greville_abscissae"       : 1,
            "parameter_resolution_for_manual_initialization"    : [ 100, 100 ],
            "max_projection_iterations"                         : 20,
            "projection_tolerance"                              : 1e-5      
        },
        "regularization_parameters":
        {
            "minimize_control_point_distance_to_surface" : true,
            "alpha"                                      : 0.001,    
            "minimize_control_point_displacement"        : false,
            "beta"                                       : 0.002
        },
        "linear_solver_name" : "SuperLU"                  
    },
    "output_parameters" : 
    {
        "output_folder"                                     : "01_Results",
        "parameter_resolution_for_output_of_surface_points" : [ 50, 50 ],
        "original_georhino_filename"                        : "tripod.georhino.txt",
        "rhino_results_filename"                            : "tripod.post.res",
        "perform_quality_evaluation"                        : true,
        "quality_evaluation_parameters" : 
        {
            "projection_parameters": 
            {
                "projection_strategy"                               : "single_search_tree",   
                "search_radius_with_multiple_trees"                 : 1.0,  
                "automatic_initialization_using_greville_abscissae" : true,
                "refinement_iterations_of_greville_abscissae"       : 1,
                "parameter_resolution_for_manual_initialization"    : [ 100, 100 ],
                "max_projection_iterations"                         : 20,
                "projection_tolerance"                              : 1e-5      
            },
            "compile_results_in_vtk_file" : true,
            "vtk_filename"                : "quality_evaluation_of_surface_reconstruction"        
        }    
    }
}""")

# ======================================================================================================================================
# Reconstruction
# ======================================================================================================================================    

print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

from cad_reconstruction_utility import CADReconstrutionUtilities

# Measure time
start_time = time.time()

# Initialize Reconstruction
CADReconstructionUtility = CADReconstrutionUtilities( reconstruction_parameters )
CADReconstructionUtility.Initialize()

# Some output before reconstruction
CADReconstructionUtility.OutputCADSurfacePoints( "surface_points_of_cad_geometry.txt" )

# Perform reconstruction
CADReconstructionUtility.PerformReconstruction()

# Some output
CADReconstructionUtility.OutputFEData()
CADReconstructionUtility.OutputCADSurfacePoints( "surface_points_of_updated_cad_geometry.txt" )
# CADReconstructionUtility.OutputGaussPointsOfFEMesh( "gauss_points_of_fe_mesh.txt" )

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")