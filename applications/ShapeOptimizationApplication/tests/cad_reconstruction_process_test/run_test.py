# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape

# Import ANurbs
import ANurbs as an

# Additional imports
import time, os, csv, shutil, sys
from contextlib import contextmanager
from cad_reconstruction_mapper import CADMapper
import numpy as np

# =======================================================================================================
# Parameters
# =======================================================================================================
specific_tests_to_run = [] # leave empty to run all tests
write_test_results_to_csv = False
suppress_output = True
delete_result_folders = True

# =======================================================================================================
# Helper functions
# =======================================================================================================
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        if suppress_output:
            sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def PerformMapping(cad_model, fe_model, parameters):
    cad_mapper = CADMapper(fe_model, cad_model, parameters)
    cad_mapper.RunMappingProcess()

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

    if write_test_results_to_csv:
        with open("pole_coords_test_"+ test_case_string +".csv", 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_NONE)
            for coords in pole_coordinates:
                csv_writer.writerow(coords)

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
# Test: Distance minimization without enforcement and alpha regularization (including rhs testing)
# =======================================================================================================
test_number = 1
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
            with open("parameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            parameters["solution"]["test_solution"].SetBool(True)
            parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

            # Mapping
            cad_model = an.Model()
            PerformMapping(cad_model, fe_model, parameters)

            # Actual test
            pole_coordinates = ExtractResults(cad_model)
            TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: Distance minimization with alpha regularization and constraint enforcement
# =======================================================================================================
test_number = 2
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["regularization"]["alpha"].SetDouble(0.1)
        parameters["conditions"]["edges"]["fe_based"]["apply_enforcement_conditions"].SetBool(True)
        parameters["conditions"]["edges"]["fe_based"]["apply_corner_enforcement_conditions"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: Distance minimization with integral approach
# =======================================================================================================
test_number = 3
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["conditions"]["general"]["apply_integral_method"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: KL shell behavior
# =======================================================================================================
test_number = 4
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: Nonlinear solution iteration (using shell that shall enforce a quasi rigid hole)
# =======================================================================================================
test_number = 5
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
            with open("parameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            parameters["solution"]["iterations"].SetInt(3)
            parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].SetBool(True)
            parameters["conditions"]["faces"]["mechanical"]["exclusive_face_list"][0].SetString("Rhino<dba9a5ba-6b33-4a1b-b46e-3916bd68d053>.BrepFace<1>.Hole")
            parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

            # Mapping
            cad_model = an.Model()
            PerformMapping(cad_model, fe_model, parameters)

            # Actual test
            pole_coordinates = ExtractResults(cad_model)
            TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: Rigid motions
# =======================================================================================================
test_number = 6
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["conditions"]["faces"]["rigid"]["apply_rigid_conditions"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: basic coupling conditions
# =======================================================================================================
test_number = 7
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["conditions"]["edges"]["coupling"]["apply_coupling_conditions"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: A prioi refinement
# =======================================================================================================
test_number = 8
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["refinement"]["a_priori"]["apply_a_priori_refinement"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Test: A posteriori refinement
# =======================================================================================================
test_number = 9
if test_number in specific_tests_to_run or len(specific_tests_to_run) == 0 :
    print("\n> Starting test "+str(test_number)+"...")

    with suppress_stdout():
        with open("parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["refinement"]["a_posteriori"]["apply_a_posteriori_refinement"].SetBool(True)
        parameters["output"]["results_directory"].SetString("Results_Test_"+str(test_number))

        # Mapping
        cad_model = an.Model()
        PerformMapping(cad_model, fe_model, parameters)

        # Actual test
        pole_coordinates = ExtractResults(cad_model)
        TestResults(pole_coordinates, str(test_number))

# =======================================================================================================
# Delete result files
# =======================================================================================================
os.remove("cad_reconstruction_process_test.post.lst")
os.remove("fe_parametrization_backup.json")

if delete_result_folders:
    for number in range(test_number+1):
        if os.path.exists("Results_Test_"+str(number)):
            shutil.rmtree("Results_Test_"+str(number))

time_for_complete_test = time.time() - start_time
print("\n> Test Successfully finished in " + str(round(time_for_complete_test,2)) + " s!\n")
# =======================================================================================================