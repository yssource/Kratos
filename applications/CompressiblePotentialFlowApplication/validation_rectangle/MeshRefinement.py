from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
#from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.CompressiblePotentialFlowApplication import *
import potential_flow_solver
import os
import loads_output_square
import shutil
#import compute_lift_process_new
import subprocess
from PyPDF2 import PdfFileReader, PdfFileMerger
from math import *
######################################################################################
######################################################################################
######################################################################################
Number_Of_Refinements = TBD

Initial_Number_Of_Segments = TBD
Refinement_Factor = TBD

work_dir = '/home/inigo/simulations/naca0012/07_salome/06_Rectangle/'
input_mdpa_path = work_dir + 'mdpas/'
output_gid_path = '/media/inigo/10740FB2740F9A1C/Outputs/06_Rectangle/'
latex_output = open(work_dir + '/plots/latex_output.txt', 'w')
latex_output.flush()

energy_h_results_file_name = work_dir + 'plots/relative_error_energy_norm/data/energy/energy_h_results.dat'
energy_n_results_file_name = work_dir + 'plots/relative_error_energy_norm/data/energy/energy_n_results.dat'
energy_variant_h_results_file_name = work_dir + 'plots/relative_error_energy_norm/data/energy/energy_variant_h_results.dat'
energy_variant_n_results_file_name = work_dir + 'plots/relative_error_energy_norm/data/energy/energy_variant_n_results.dat'
energy_results_directory_name = work_dir + 'plots/relative_error_energy_norm/data/energy'

condition_results_file_name = work_dir + 'plots/condition_number/data/condition/condition_results.dat'
condition_results_directory_name = work_dir + 'plots/condition_number/data/condition'

loads_output_square.write_header_all_cases(work_dir)

merger_global_far_field_x = PdfFileMerger()
merger_global_far_field_y = PdfFileMerger()

case = 0

merger_local_far_field_x = PdfFileMerger()
merger_local_far_field_y = PdfFileMerger()

with open(energy_h_results_file_name,'w') as energy_file:
    energy_file.flush()

with open(energy_n_results_file_name,'w') as energy_file:
    energy_file.flush()

with open(energy_variant_h_results_file_name,'w') as energy_file:
    energy_file.flush()

with open(energy_variant_n_results_file_name,'w') as energy_file:
    energy_file.flush()

with open(condition_results_file_name, 'w') as condition_file:
    condition_file.flush()

AOA = 5.0
MeshSize = 1.0
mesh_refinement_file_name = work_dir + 'plots/results/mesh_refinement'
energy_data_directory_name = 'data/energy_results'
condition_data_directory_name = 'data/condition_results'
loads_output_square.write_header(work_dir)

far_field_data_directory_start = work_dir + 'plots/far_field/data/far_field'
os.mkdir(far_field_data_directory_start)

os.mkdir(output_gid_path + 'gid_results')

energy_reference = 1.0
potential_energy_reference = 1.0
counter = 0.0

NumberOfSegments_tmp = Initial_Number_Of_Segments

for i in range(Number_Of_Refinements):
    NumberOfSegments = int(NumberOfSegments_tmp)
    MeshSize = 100.0/NumberOfSegments
    print("\n\tCase ", case, "\n")
    loads_output_square.write_case(case, NumberOfSegments, MeshSize, work_dir)

    far_field_results_directory_name = work_dir + 'plots/far_field/data/0_original'
    far_field_data_directory_name = far_field_data_directory_start + '/Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments) 

    ## Parse the ProjectParameters
    #parameter_file = open("ProjectParameters_compressibility.json",'r')
    parameter_file = open("ProjectParameters.json",'r')
    ProjectParameters = Parameters( parameter_file.read())

    ## Get echo level and parallel type
    verbosity = ProjectParameters["problem_data"]["echo_level"].GetInt()
    parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

    ## Fluid model part definition
    main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
    main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

    ###TODO replace this "model" for real one once available
    Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

    #Set Mesh input_filename
    mdpa_file_name = input_mdpa_path + 'naca0012_Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments)

    ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mdpa_file_name)
    ProjectParameters["boundary_conditions_process_list"][1]["Parameters"]["meshsize"].SetDouble(MeshSize)
    ProjectParameters["boundary_conditions_process_list"][1]["Parameters"]["energy_reference"].SetDouble(energy_reference)
    ProjectParameters["boundary_conditions_process_list"][1]["Parameters"]["potential_energy_reference"].SetDouble(potential_energy_reference)

    ## Solver construction    
    solver = potential_flow_solver.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

    solver.AddVariables()

    ## Read the model - note that SetBufferSize is done here
    solver.ImportModelPart()

    ## Add AddDofs
    solver.AddDofs()

    ## Set output name
    #problem_name = output_gid_path + ProjectParameters["problem_data"]["problem_name"].GetString()+ "Mesh" + str(case)
    problem_name = output_gid_path + 'gid_results/' + ProjectParameters["problem_data"]["problem_name"].GetString()+ '_Case_' + str(
        case)+ '_NumberOfSegments_' + str(NumberOfSegments)

    ## Initialize GiD  I/O
    if (parallel_type == "OpenMP"):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                    problem_name ,
                                    ProjectParameters["output_configuration"])

    gid_output.ExecuteInitialize()

    ##TODO: replace MODEL for the Kratos one ASAP
    ## Get the list of the skin submodel parts in the object Model
    for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
        skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
        Model.update({skin_part_name: main_model_part.GetSubModelPart(skin_part_name)})

    ## Get the list of the no-skin submodel parts in the object Model (results processes and no-skin conditions)
    for i in range(ProjectParameters["solver_settings"]["no_skin_parts"].size()):
        no_skin_part_name = ProjectParameters["solver_settings"]["no_skin_parts"][i].GetString()
        Model.update({no_skin_part_name: main_model_part.GetSubModelPart(no_skin_part_name)})

    ## Print model_part and properties
    if(verbosity > 1):
        print("")
        print(main_model_part)
        for properties in main_model_part.Properties:
            print(properties)

    ## Processes construction
    import process_factory
    # "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
    # Note 1: gravity is firstly constructed. Outlet process might need its information.
    # Note 2: conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
    #ProjectParameters["boundary_conditions_process_list"]["model_import_settings"]["input_filename"].SetString(mdpa_file_name)
    #ProjectParameters["boundary_conditions_process_list"][2]["Parameters"]["mesh_refinement_file_name"].SetString(mesh_refinement_file_name)
    list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

    if(verbosity > 1):
        for process in list_of_processes:
            print(process)

    ## Processes initialization
    for process in list_of_processes:
        process.ExecuteInitialize()

    ## Solver initialization
    solver.Initialize()

    #TODO: think if there is a better way to do this
    fluid_model_part = solver.GetComputingModelPart()

    ## Stepping and time settings
    # Dt = ProjectParameters["problem_data"]["time_step"].GetDouble()
    start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
    end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

    time = start_time
    step = 0
    out = 0.0

    gid_output.ExecuteBeforeSolutionLoop()

    for process in list_of_processes:
        process.ExecuteBeforeSolutionLoop()

    ## Writing the full ProjectParameters file before solving
    if ((parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0)) and (verbosity > 0):
        f = open("ProjectParametersOutput.json", 'w')
        f.write(ProjectParameters.PrettyPrintJsonString())
        f.close()

    Dt = 0.01
    step += 1
    time = time + Dt
    main_model_part.CloneTimeStep(time)

    if (parallel_type == "OpenMP") or (KratosMPI.mpi.rank == 0):
        print("")
        print("STEP = ", step)
        print("TIME = ", time)

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()

    solver.Solve()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    gid_output.ExecuteFinalizeSolutionStep()

    #TODO: decide if it shall be done only when output is processed or not
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    if gid_output.IsOutputStep():
        gid_output.PrintOutput()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    out = out + Dt

    gid_output.ExecuteFinalize()

    #output far field
    loads_output_square.write_figures_far_field(far_field_data_directory_name, AOA, case, MeshSize,  NumberOfSegments, work_dir)
    shutil.copytree(far_field_results_directory_name, far_field_data_directory_name)

    latex_far_field = subprocess.Popen(['pdflatex', '-interaction=batchmode',work_dir + 'plots/far_field/main_far_field_x.tex'], stdout=latex_output)
    latex_far_field.communicate()

    latex_far_field = subprocess.Popen(['pdflatex', '-interaction=batchmode',work_dir + 'plots/far_field/main_far_field_y.tex'], stdout=latex_output)
    latex_far_field.communicate()

    far_field_x_file_name = work_dir + 'plots/far_field/plots/velocity_x_Case_' + str(case) + '_AOA_' + str(
            AOA) + '_Far_Field_Mesh_Size_' + str(NumberOfSegments) + '_Airfoil_Mesh_Size_' + str(MeshSize) + '.pdf'
    shutil.copyfile('main_far_field_x.pdf',far_field_x_file_name)
    merger_local_far_field_x.append(PdfFileReader(far_field_x_file_name), 'case_' + str(case))
    merger_global_far_field_x.append(PdfFileReader(far_field_x_file_name), 'case_' + str(case))

    far_field_y_file_name = work_dir + 'plots/far_field/plots/velocity_y_Case_' + str(case) + '_AOA_' + str(
            AOA) + '_Far_Field_Mesh_Size_' + str(NumberOfSegments) + '_Airfoil_Mesh_Size_' + str(MeshSize) + '.pdf'
    shutil.copyfile('main_far_field_y.pdf',far_field_y_file_name)
    merger_local_far_field_y.append(PdfFileReader(far_field_y_file_name), 'case_' + str(case))
    merger_global_far_field_y.append(PdfFileReader(far_field_y_file_name), 'case_' + str(case))

    if(counter < 1):
        energy_reference = main_model_part.ProcessInfo[ENERGY_NORM_REFERENCE]
        potential_energy_reference = main_model_part.ProcessInfo[POTENTIAL_ENERGY_REFERENCE]

    NumberOfSegments_tmp *= Refinement_Factor

    case +=1
    counter +=1.0
    
for process in list_of_processes:
    process.ExecuteFinalize()

os.rename(work_dir + "mesh_refinement_loads.dat", mesh_refinement_file_name)

loads_output_square.write_figures_condition(condition_data_directory_name, AOA, work_dir)
loads_output_square.write_figures_energy(energy_data_directory_name, AOA, work_dir)

shutil.copytree(energy_results_directory_name, work_dir + 'plots/relative_error_energy_norm/' + energy_data_directory_name)
os.remove(energy_h_results_file_name)
os.remove(energy_n_results_file_name)
os.remove(energy_variant_h_results_file_name)
os.remove(energy_variant_n_results_file_name)

shutil.copytree(condition_results_directory_name, work_dir + 'plots/condition_number/' + condition_data_directory_name)
os.remove(condition_results_file_name)

far_field_x_final_file_name = work_dir + 'plots/far_field/far_field_x_AOA_' + str(AOA) + '.pdf'
merger_local_far_field_x.write(far_field_x_final_file_name)

far_field_y_final_file_name = work_dir + 'plots/far_field/far_field_y_AOA_' + str(AOA) + '.pdf'
merger_local_far_field_y.write(far_field_y_final_file_name)


far_field_x_global_file_name = work_dir + 'plots/far_field/far_field_x_all.pdf'
merger_global_far_field_x.write(far_field_x_global_file_name)

far_field_y_global_file_name = work_dir + 'plots/far_field/far_field_y_all.pdf'
merger_global_far_field_y.write(far_field_y_global_file_name)
