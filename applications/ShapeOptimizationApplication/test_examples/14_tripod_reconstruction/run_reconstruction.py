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

# Input parameters
fem_filename = "tripod"
cad_geometry_filename = "tripod_geometry.json" 
cad_integration_data_filename = "tripod_integration_data.json"

# Output parameters
results_output_settings = Parameters("""
{
    "results_output_folder"                             : "01_Results",
    "parameter_resolution_for_output_of_surface_points" : [ 50, 50 ],
    "original_georhino_filename"                        : "tripod.georhino.txt",
    "rhino_results_filename"                            : "tripod.post.res"
}""")

# ======================================================================================================================================
# Reconstruction
# ======================================================================================================================================    

print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

import cad_reconstruction_utility

# Measure time
start_time = time.time()

# Initialize Reconstruction
CADReconstructionUtility = cad_reconstruction_utility.CADReconstrutionUtilities( fem_filename, cad_geometry_filename, cad_integration_data_filename, results_output_settings )
CADReconstructionUtility.Initialize()

# Set Boundary Conditions
CADReconstructionUtility.SetDisplacementCouplingOnAllCouplingPoints()
CADReconstructionUtility.SetRotationCouplingOnAllCouplingPoints()

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