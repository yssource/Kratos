# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape

# Import ANurbs
import ANurbs as an

# Additional imports
import time, os, csv, shutil, sys
from contextlib import contextmanager
from cad_reconstruction import CADMapper
import numpy as np

# =======================================================================================================
# Helper functions
# =======================================================================================================
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def PerformMapping(cad_model, fe_model):
    cad_mapper = CADMapper(fe_model, cad_model, parameters)
    cad_mapper.ReadModelData()
    cad_mapper.Initialize()
    cad_mapper.Map()
    cad_mapper.Finalize()

def ExtractResults(cad_model):
    pole_coordinates = []
    for surface_i in cad_model.GetByType('SurfaceGeometry3D'):
        surface_geometry = surface_i.Data()
        surface_geometry_key = surface_i.Key()

        for r in range(surface_geometry.NbPolesU()):
            for s in range(surface_geometry.NbPolesV()):
                pole_coordinates.append(surface_geometry.Pole(r,s))
    return pole_coordinates

def TestResults(pole_coordinates, test_case_string):
    error_tolerance = 10

    # with open("pole_coords_test_"+ test_case_string +".csv", 'w', newline='') as csvfile:
    #     csv_writer = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_NONE)
    #     for coords in pole_coordinates:
    #         csv_writer.writerow(coords)

    with open("pole_coords_test_"+ test_case_string +".csv", 'r') as f:
        reader = csv.reader(f)
        reference_coordinates = []
        for row in reader:
            reference_coordinates.append( np.array([ float(value) for value in row ]) )
    np.testing.assert_array_almost_equal(np.array(pole_coordinates), np.array(reference_coordinates), decimal=error_tolerance)

# =======================================================================================================
# Preprocessing
# =======================================================================================================
fe_model = KratosMultiphysics.Model()
start_time = time.time()

# =======================================================================================================
# Test 1: Distance minimization without enforcement and alpha regularization (including several iterations and rhs testing)
# =======================================================================================================
print("\n> Starting test 1...")
start_time_2 = time.time()

with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["solution"]["iterations"].SetInt(3)
    parameters["solution"]["test_solution"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_1")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "01")

time_for_first_test = time.time()-start_time_2

# =======================================================================================================
# Test 2: Distance minimization with alpha regularization and constraint enforcement
# =======================================================================================================
print("\n> Starting test 2...")
with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["regularization"]["alpha"].SetDouble(0.1)
    parameters["conditions"]["edges"]["fe_based"]["apply_enforcement_conditions"].SetBool(True)
    parameters["conditions"]["edges"]["fe_based"]["apply_corner_enforcement_conditions"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_2")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "02")

# =======================================================================================================
# Test 3: Distance minimization with integral approach
# =======================================================================================================
print("\n> Starting test 3...")
with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["conditions"]["general"]["apply_integral_method"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_3")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "03")

# =======================================================================================================
# Test 4: Curvature condition
# =======================================================================================================
print("\n> Starting test 4...")
with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["conditions"]["faces"]["curvature"]["apply_curvature_minimization"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_4")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "04")

# =======================================================================================================
# Test 5: KL shell behavior
# =======================================================================================================
print("\n> Starting test 5...")
with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_5")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "05")

# =======================================================================================================
# Test 6: Rigid motions
# =======================================================================================================
print("\n> Starting test 6...")
with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["conditions"]["faces"]["rigid"]["apply_rigid_conditions"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_6")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "06")

# =======================================================================================================
# Test 7: basic coupling conditions
# =======================================================================================================
print("\n> Starting test 7...")
with suppress_stdout():

    with open("parameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters["conditions"]["edges"]["coupling"]["apply_coupling_conditions"].SetBool(True)
    parameters["output"]["results_directory"].SetString("Results_Test_7")

    # Mapping
    cad_model = an.Model()
    PerformMapping(cad_model, fe_model)

    # Actual test
    pole_coordinates = ExtractResults(cad_model)
    TestResults(pole_coordinates, "07")

# =======================================================================================================
# Test time
# =======================================================================================================
print("\n> Starting time...")
time_for_complete_test = time.time() - start_time
relative_time_ratio = time_for_complete_test / time_for_first_test

reference_ratio = 11.2

if (relative_time_ratio - reference_ratio) / relative_time_ratio > 0.10:
    raise RuntimeError("Test took unexpectedly long!")
elif (relative_time_ratio - reference_ratio) / relative_time_ratio < -0.10:
    raise RuntimeError("Test took unexpectedly short!")

# =======================================================================================================
# Delete result files
# =======================================================================================================
os.remove("cad_reconstruction_process_test.post.lst")
os.remove("plate.time")
shutil.rmtree("Results_Test_1")
shutil.rmtree("Results_Test_2")
shutil.rmtree("Results_Test_3")
shutil.rmtree("Results_Test_4")
shutil.rmtree("Results_Test_5")
shutil.rmtree("Results_Test_6")
shutil.rmtree("Results_Test_7")

print("\n> Test Successfully finished in " + str(round(time_for_complete_test,2)) + " s!\n")
# =======================================================================================================