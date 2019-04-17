# Import Kratos core and apps
import KratosMultiphysics as km
import KratosMultiphysics.ShapeOptimizationApplication as kso
import KratosMultiphysics.StructuralMechanicsApplication as kcsm
import KratosMultiphysics.MappingApplication as kma
import structural_response_function_factory

# Additional imports
from analyzer_base import AnalyzerBaseClass
from gid_output_process import GiDOutputProcess
import time, os, shutil, csv, math
from analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication import custom_variable_utilities as cvu

# =================================================================================================================================================================
# Parameters
# =================================================================================================================================================================
# verification parameters
fd_start_level = 0
fd_end_level = 10
list_of_move_nodes = [2114]

csm_point_load_filename = "csm_point_load_baseline_non_dimensional.txt"

stress_response_settings = km.Parameters("""
{
        "response_type"        : "adjoint_max_stress",
        "gradient_mode"        : "semi_analytic",
        "step_size"            : 1e-11,
        "critical_part_name"   : "stress_partition",
        "stress_type"          : "VON_MISES_STRESS",
        "stress_treatment"     : "mean",
        "echo_level"           : 1,
        "primal_settings"      : "parameters_analysis.json",
        "adjoint_settings"     : "auto",
        "sensitivity_settings" : {
            "sensitivity_model_part_name"     : "Parts_structure",
            "nodal_sensitivity_variables"     : ["SHAPE"],
            "element_sensitivity_variables"   : [],
            "condition_sensitivity_variables" : [],
            "build_mode": "static"
        }
}""")

with open("parameters_optimization.json",'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())

# =================================================================================================================================================================
# Helper functions
# =================================================================================================================================================================
# -------------------------------------------------------------------------------------------------------------------------------------------------------------
def OutputOPTAsGid( output_mdpa, output_filename ):
    output_parameters = km.Parameters("""
    {
        "result_file_configuration" : {
            "gidpost_flags"       : {
                "GiDPostMode"           : "GiD_PostBinary",
                "WriteDeformedMeshFlag" : "WriteDeformed",
                "WriteConditionsFlag"   : "WriteConditions",
                "MultiFileFlag"         : "SingleFile"
            },
            "file_label"          : "step",
            "output_control_type" : "step",
            "output_frequency"    : 1,
            "body_output"         : true,
            "node_output"         : false,
            "skin_output"         : false,
            "nodal_results"       : ["CONTROL_POINT_UPDATE","SHAPE_UPDATE"],
            "gauss_point_results" : []
        },
        "point_data_configuration"  : []
    }""")
    gid_output_original = GiDOutputProcess( output_mdpa, output_filename, output_parameters )
    gid_output_original.ExecuteInitialize()
    gid_output_original.ExecuteBeforeSolutionLoop()
    gid_output_original.ExecuteInitializeSolutionStep()
    gid_output_original.PrintOutput()
    gid_output_original.ExecuteFinalizeSolutionStep()
    gid_output_original.ExecuteFinalize()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
def OutputOPTAsGidWithGradients( output_mdpa, output_filename ):
    output_parameters = km.Parameters("""
    {
        "result_file_configuration" : {
            "gidpost_flags"       : {
                "GiDPostMode"           : "GiD_PostBinary",
                "WriteDeformedMeshFlag" : "WriteDeformed",
                "WriteConditionsFlag"   : "WriteConditions",
                "MultiFileFlag"         : "SingleFile"
            },
            "file_label"          : "step",
            "output_control_type" : "step",
            "output_frequency"    : 1,
            "body_output"         : true,
            "node_output"         : false,
            "skin_output"         : false,
            "nodal_results"       : ["CONTROL_POINT_UPDATE","SHAPE_UPDATE","DF1DX","DF1DX_MAPPED"],
            "gauss_point_results" : []
        },
        "point_data_configuration"  : []
    }""")
    gid_output_original = GiDOutputProcess( output_mdpa, output_filename, output_parameters )
    gid_output_original.ExecuteInitialize()
    gid_output_original.ExecuteBeforeSolutionLoop()
    gid_output_original.ExecuteInitializeSolutionStep()
    gid_output_original.PrintOutput()
    gid_output_original.ExecuteFinalizeSolutionStep()
    gid_output_original.ExecuteFinalize()

# =================================================================================================================================================================
# Define analyzer
# =================================================================================================================================================================
class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, calculate_gradient=False):
        # obtain new csm mesh
        new_mesh = {node.Id: [node.X, node.Y, node.Z] for node in current_design.Nodes}

        # Perform csm
        print("\n> Starting __RunCSM in analyzer...")
        start_time = time.time()

        analysis_model = km.Model()

        csm = structural_response_function_factory.CreateResponseFunction(stress_response_settings["response_type"].GetString(), stress_response_settings, analysis_model)

        csm_primal_part = csm.primal_model_part
        csm_primal_part.AddNodalSolutionStepVariable(kso.MESH_CHANGE)
        csm_primal_part.AddNodalSolutionStepVariable(kso.TRACTION)

        # Initialize and set time to current optimization iteration
        csm.primal_analysis.project_parameters["problem_data"]["start_time"].SetDouble(optimization_iteration-1)
        csm.primal_analysis.project_parameters["problem_data"]["end_time"].SetDouble(optimization_iteration)

        csm.Initialize()

        csm_primal_part.ProcessInfo.SetValue(km.STEP, optimization_iteration-1)

        # apply mesh motion
        if new_mesh != {}:
            for node in csm_primal_part.Nodes:
                X_new = new_mesh[node.Id]
                mesh_change = [X_new[0]-node.X, X_new[1]-node.Y, X_new[2]-node.Z]
                node.SetSolutionStepValue(kso.MESH_CHANGE, mesh_change)
                node.X = X_new[0]
                node.Y = X_new[1]
                node.Z = X_new[2]
                node.X0 = X_new[0]
                node.Y0 = X_new[1]
                node.Z0 = X_new[2]

        # compute primal field

        # Apply forces
        # Point load is read from file which contains results from a cfd simulation of the undeformed design up to a residual of -15
        with open(csm_point_load_filename, 'r') as infile:
            csv_reader = csv.reader(infile, delimiter=',')
            lines = list(csv_reader)
        point_loads = {int(line[0]): [float(line[1]), float(line[2]), float(line[3])] for line in lines}
        for node in csm_primal_part.Nodes:
            node.SetSolutionStepValue(kcsm.POINT_LOAD, point_loads[node.Id])

        # Compute primals
        csm.InitializeSolutionStep()
        csm.CalculateValue()
        if calculate_gradient:
            csm.CalculateGradient()
        csm.FinalizeSolutionStep()
        csm.Finalize()

        print("> Finished __RunCSM in" ,round( time.time()-start_time, 3 ), " s.")

        # Return
        if calculate_gradient:
            return csm.GetValue(), csm.GetShapeGradient()
        else:
            return csm.GetValue()

# =================================================================================================================================================================
# Perform verification
# =================================================================================================================================================================
for move_node_id in list_of_move_nodes:
    # Compute reference
    print("\n\n##########################################################################")
    print(">> Starting to compute reference of node "+str(move_node_id))
    print("##########################################################################")
    start_time = time.time()

    # Create case folder
    case_dir = "tmp_calculation_reference"
    if os.path.exists(case_dir):
        os.system("rm -rf "+case_dir+"/")
    os.system("mkdir "+case_dir)

    # Copy all necessary files to case folder
    shutil.copy("CSM_model.mdpa", case_dir+"/")
    shutil.copy("OPT_model.mdpa", case_dir+"/")
    shutil.copy("parameters_analysis.json", case_dir+"/")
    shutil.copy("parameters_material.json", case_dir+"/")
    shutil.copy("parameters_optimization.json", case_dir+"/")
    shutil.copy(csm_point_load_filename, case_dir+"/")

    # Change to case directory
    top_dir = os.getcwd()
    os.chdir(case_dir)

    # Create analyzer
    analyzer = CustomAnalyzer()

    # Read optimization model part
    model = km.Model()
    opt_part = model.CreateModelPart("OPT_model")
    opt_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)
    opt_part.AddNodalSolutionStepVariable(kso.CONTROL_POINT_UPDATE)
    opt_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)
    opt_part.AddNodalSolutionStepVariable(kso.DF1DX)
    opt_part.AddNodalSolutionStepVariable(kso.DF1DX_MAPPED)

    model_part_io = km.ModelPartIO("OPT_model")
    model_part_io.ReadModelPart(opt_part)

    # Run analysis
    analyzer.InitializeBeforeOptimizationLoop()
    csm_response_value_reference, csm_response_gradient_reference = analyzer.AnalyzeDesignAndReportToCommunicator(opt_part, 1, calculate_gradient=True)

    # Store gradient on node
    for node in opt_part.Nodes:
        node.SetSolutionStepValue(kso.DF1DX, csm_response_gradient_reference[node.Id])

    # Map gradient
    vm_mapper = kso.MapperVertexMorphingMatrixFree(opt_part, opt_part, parameters["optimization_settings"]["design_variables"]["filter"])
    vm_mapper.InverseMap(kso.DF1DX, kso.DF1DX_MAPPED)

    # Some output
    OutputOPTAsGidWithGradients(opt_part, "opt_part")

    # Save values in file
    with open("csm_value.txt",'w') as results_file:
        results_file.write(str(csm_response_value_reference))

    # Change back to top dir and save results
    os.chdir(top_dir)

    with open("fd_results_wo_fluid_move_node_"+str(move_node_id)+".txt",'w') as results_file:
        results_file.write("csm_reference_value = "+str(csm_response_value_reference)+"\n")
        results_file.write("csm_reference_gradient = "+str(opt_part.Nodes[move_node_id].GetSolutionStepValue(kso.DF1DX_MAPPED))+"\n")
        results_file.write("---\n")
        results_file.write("step_size,csm_gradient_x,csm_gradient_y,time\n")

    print("\n> Finished computing reference value in" ,round( time.time()-start_time, 3 ), " s.")
    print("> csm reference response value =",csm_response_value_reference)
    print("> csm reference response gradient [move_node_id] =",opt_part.Nodes[move_node_id].GetSolutionStepValue(kso.DF1DX_MAPPED))
    print("\n")


    # Compute pertubation
    for itr in range(fd_start_level,fd_end_level+1):

        # Compute pertubation
        current_delta = 1*10**(-itr)

        for dim_itr in range(2):

            print("\n\n##########################################################################")
            print(">> Starting finite difference level ",itr, " dimension ",dim_itr+1, "for node ", str(move_node_id))
            print("##########################################################################")

            # Create case folder
            case_dir = "tmp_calculation_"+str(itr)+"_"+str(dim_itr+1)
            if os.path.exists(case_dir):
                os.system("rm -rf "+case_dir+"/")
            os.system("mkdir "+case_dir)

            # Copy all necessary files to case folder
            shutil.copy("CSM_model.mdpa", case_dir+"/")
            shutil.copy("OPT_model.mdpa", case_dir+"/")
            shutil.copy("parameters_analysis.json", case_dir+"/")
            shutil.copy("parameters_material.json", case_dir+"/")
            shutil.copy("parameters_optimization.json", case_dir+"/")
            shutil.copy(csm_point_load_filename, case_dir+"/")

            # Change to case directory
            top_dir = os.getcwd()
            os.chdir(case_dir)

            # Create analyzer
            analyzer = CustomAnalyzer()

            # Read optimization model part
            model = km.Model()
            opt_part = model.CreateModelPart("OPT_model")
            opt_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)
            opt_part.AddNodalSolutionStepVariable(kso.CONTROL_POINT_UPDATE)
            opt_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)

            model_part_io = km.ModelPartIO("OPT_model")
            model_part_io.ReadModelPart(opt_part)

            # Apply control point update
            if dim_itr == 0:
                opt_part.Nodes[move_node_id].SetSolutionStepValue(kso.CONTROL_POINT_UPDATE, [current_delta,0,0])
            elif dim_itr == 1:
                opt_part.Nodes[move_node_id].SetSolutionStepValue(kso.CONTROL_POINT_UPDATE, [0,current_delta,0])

            # Compute corresponding shape update
            vm_mapper = kso.MapperVertexMorphingMatrixFree(opt_part, opt_part, parameters["optimization_settings"]["design_variables"]["filter"])
            vm_mapper.Map(kso.CONTROL_POINT_UPDATE, kso.SHAPE_UPDATE)

            # Update mesh
            kso.MeshControllerUtilities(opt_part).UpdateMeshAccordingInputVariable(kso.SHAPE_UPDATE)
            kso.MeshControllerUtilities(opt_part).SetReferenceMeshToMesh()

            # Run analysis
            analyzer.InitializeBeforeOptimizationLoop()
            csm_response_value = analyzer.AnalyzeDesignAndReportToCommunicator(opt_part, 1, calculate_gradient=False)

            # Evaluate fd gradients
            if dim_itr == 0:
                csm_grad_x = (csm_response_value - csm_response_value_reference) / current_delta
            elif dim_itr == 1:
                csm_grad_y = (csm_response_value - csm_response_value_reference) / current_delta

            print("> current csm response value =",csm_response_value)
            print("\n")

            # Change back to top dir
            os.chdir(top_dir)

        # Write results
        with open("fd_results_wo_fluid_move_node_"+str(move_node_id)+".txt",'a') as results_file:
            results_file.write(str(current_delta)+","+str(csm_grad_x)+","+str(csm_grad_y)+","+str(time.ctime())+"\n")

    print("\n\n##########################################################################")
    print("> Finished verifiation process in" ,round( time.time()-start_time, 3 ), " s.")
    print("##########################################################################")
