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
fd_start_level = 4
fd_end_level = 15
list_of_test_nodes = [3679]

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
    def AnalyzeDesignAndReportToCommunicator(self, optimization_iteration, calculate_gradient=False):

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
start_time = time.time()

for test_node_id in list_of_test_nodes:
    with open("fd_semianalytic_csm_gradient_test_node_"+str(test_node_id)+".txt",'w') as results_file:
        results_file.write("step_size,csm_gradient_x,csm_gradient_y,time\n")

    # Compute pertubation
    for itr in range(fd_start_level,fd_end_level+1):

        # Compute pertubation
        current_delta = 1*10**(-itr)

        print("\n\n##########################################################################")
        print(">> Starting finite difference level ",itr, "for node ", str(test_node_id))
        print("##########################################################################")

        # Create case folder
        case_dir = "tmp_calculation"
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

        # Apply pertubation
        stress_response_settings["step_size"].SetDouble(current_delta)

        # Run analysis
        analyzer.InitializeBeforeOptimizationLoop()
        csm_response_value, csm_gradient = analyzer.AnalyzeDesignAndReportToCommunicator(1, calculate_gradient=True)

        # Store gradient on node
        for node in opt_part.Nodes:
            node.SetSolutionStepValue(kso.DF1DX, csm_gradient[node.Id])

        # Map gradient
        vm_mapper = kso.MapperVertexMorphingMatrixFree(opt_part, opt_part, parameters["optimization_settings"]["design_variables"]["filter"])
        vm_mapper.InverseMap(kso.DF1DX, kso.DF1DX_MAPPED)

        csm_grad_of_test_node = opt_part.Nodes[test_node_id].GetSolutionStepValue(kso.DF1DX_MAPPED)

        print("> current csm response value =",csm_response_value)
        print("> current csm response gradient ="+str(csm_grad_of_test_node)+"\n")

        print("\n")

        # Change back to top dir
        os.chdir(top_dir)

        # Write results
        with open("fd_semianalytic_csm_gradient_test_node_"+str(test_node_id)+".txt",'a') as results_file:
            results_file.write(str(current_delta)+","+str(csm_grad_of_test_node[0])+","+str(csm_grad_of_test_node[1])+","+str(time.ctime())+"\n")

print("\n\n##########################################################################")
print("> Finished verifiation process in" ,round( time.time()-start_time, 3 ), " s.")
print("##########################################################################")
