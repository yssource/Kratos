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
move_node_id = 3679

response_function = "DRAG"

# analysis parameter
cfd_part_name = "ubend_2d"
cfd_interface_part_name = "design_surface"
cfd_wall_part_name = "WALL"

with open("parameters_optimization.json",'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())

# Change threading layer, otherwise the suprocess module used in the interface.py in SU2 will hang, (only)
# when Intel solvers are used as linear solvers for calculating the structure. This is a known issue of
# the Intel compiler 2018 and will be fixed in future releases. Compare:
# 1) https://github.com/numpy/numpy/issues/10060
# 2) https://software.intel.com/en-us/forums/intel-c-compiler/topic/758961
os.environ["MKL_THREADING_LAYER"] = "TBB"

# =================================================================================================================================================================
# Helper functions
# =================================================================================================================================================================
def OutputCFDAsGid( output_mdpa, output_filename ):
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
            "nodal_results"       : ["SHAPE_UPDATE","NORMALIZED_SURFACE_NORMAL", "DF1DX"],
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
            "nodal_results"       : ["CONTROL_POINT_UPDATE","SHAPE_UPDATE","DF1DX_MAPPED"],
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
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __init__(self):
        self.cfd_interface_part = None
        self.cfd_part = None
        self.cfd_wall_part = None

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        # Initialize cfd data
        cfd_model = km.Model()
        self.cfd_part = cfd_model.CreateModelPart(cfd_part_name)
        self.cfd_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,2)
        self.cfd_part.AddNodalSolutionStepVariable(kso.TRACTION)
        self.cfd_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)
        self.cfd_part.AddNodalSolutionStepVariable(km.NORMAL)
        self.cfd_part.AddNodalSolutionStepVariable(kso.DF1DX)
        self.cfd_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)

        model_part_io = km.ModelPartIO(cfd_part_name)
        model_part_io.ReadModelPart(self.cfd_part)

        self.cfd_interface_part = self.cfd_part.GetSubModelPart(cfd_interface_part_name)
        self.cfd_wall_part = self.cfd_part.GetSubModelPart(cfd_wall_part_name)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, calculate_gradient=False):
        # obtain new csm mesh
        new_csm_mesh = {node.Id: [node.X, node.Y, node.Z] for node in current_design.Nodes}

        # Update fluid mesh (structure is controlled by the optimization algorithm) - Note that the mapper should be based on the previos design to match the forward map in the structure
        for node in current_design.Nodes:
            shape_update = node.GetSolutionStepValue(kso.SHAPE_UPDATE)
            node.X -= shape_update[0]
            node.Y -= shape_update[1]
            node.Z -= shape_update[2]
            node.X0 -= shape_update[0]
            node.Y0 -= shape_update[1]
            node.Z0 -= shape_update[2]

        vm_mapper = kso.MapperVertexMorphingMatrixFree(current_design, self.cfd_interface_part, parameters["optimization_settings"]["design_variables"]["filter"])
        vm_mapper.Map(kso.CONTROL_POINT_UPDATE, kso.SHAPE_UPDATE)

        for node in current_design.Nodes:
            shape_update = node.GetSolutionStepValue(kso.SHAPE_UPDATE)
            node.X += shape_update[0]
            node.Y += shape_update[1]
            node.Z += shape_update[2]
            node.X0 += shape_update[0]
            node.Y0 += shape_update[1]
            node.Z0 += shape_update[2]

        kso.MeshControllerUtilities(self.cfd_interface_part).UpdateMeshAccordingInputVariable(kso.SHAPE_UPDATE)
        kso.MeshControllerUtilities(self.cfd_interface_part).SetReferenceMeshToMesh()

        # Initialize SU2 interface
        from interface_su2 import InterfaceSU2
        self.interface_su2 = InterfaceSU2(parameters["su2_interface_settings"])
        self.interface_su2.InitializeNewSU2Project()

        self.interface_su2.WriteNodesAsSU2MeshMotionFile(self.cfd_wall_part.GetNodes())

        # Run fluid
        update_mesh = True
        [cfd_response_value] = self.interface_su2.ComputeValues([response_function], update_mesh, optimization_iteration)

        if calculate_gradient:
            update_mesh = False
            [cfd_response_gradient] = self.interface_su2.ComputeGradient([response_function], update_mesh, optimization_iteration)
            cvu.WriteDictionaryDataOnNodalVariable(cfd_response_gradient, self.cfd_part, kso.DF1DX)

        # Some output
        OutputCFDAsGid(self.cfd_interface_part, "cfd_interface_part")

        # Return
        return cfd_response_value

# =================================================================================================================================================================
# Perform verification
# =================================================================================================================================================================
# Compute reference
print("\n\n##########################################################################")
print(">> Starting to compute reference")
print("##########################################################################")
start_time = time.time()

# Create case folder
case_dir = "tmp_calculation_reference"
if os.path.exists(case_dir):
    os.system("rm -rf "+case_dir+"/")
os.system("mkdir "+case_dir)

# Copy all necessary files to case folder
shutil.copy("ubend_2d.mdpa", case_dir+"/")
shutil.copy("OPT_model.mdpa", case_dir+"/")
shutil.copy("parameters_optimization.json", case_dir+"/")
shutil.copy("restart_flow.dat", case_dir+"/")
shutil.copy("ubend_coarse.cfg", case_dir+"/")
shutil.copy("ubend_2d.su2", case_dir+"/")
shutil.copy("su2_mesh_data.json", case_dir+"/")

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
opt_part.AddNodalSolutionStepVariable(kso.DF1DX_MAPPED)

model_part_io = km.ModelPartIO("OPT_model")
model_part_io.ReadModelPart(opt_part)

# Run analysis
analyzer.InitializeBeforeOptimizationLoop()
cfd_response_value_reference = analyzer.AnalyzeDesignAndReportToCommunicator(opt_part, 1, calculate_gradient=True)

# Map gradient
vm_mapper = kso.MapperVertexMorphingMatrixFree(opt_part, analyzer.cfd_part, parameters["optimization_settings"]["design_variables"]["filter"])
vm_mapper.InverseMap(kso.DF1DX, kso.DF1DX_MAPPED)

# Some output
OutputOPTAsGidWithGradients(opt_part, "opt_part")

# Save values in file
with open("cfd_value.txt",'w') as results_file:
    results_file.write(str(cfd_response_value_reference))

# Change back to top dir and save results
os.chdir(top_dir)

with open("fd_results_"+str(response_function)+"_move_node_"+str(move_node_id)+".txt",'w') as results_file:
    results_file.write("cfd_reference_value = "+str(cfd_response_value_reference)+"\n")
    results_file.write("cfd_reference_gradient = "+str(opt_part.Nodes[move_node_id].GetSolutionStepValue(kso.DF1DX_MAPPED))+"\n")
    results_file.write("---\n")
    results_file.write("step_size,cfd_gradient_x,cfd_gradient_y,time\n")

print("\n> Finished computing reference value in" ,round( time.time()-start_time, 3 ), " s.")
print("> cfd reference response value =",cfd_response_value_reference)
print("\n")

# Compute pertubation
for itr in range(fd_start_level,fd_end_level+1):

    # Compute pertubation
    current_delta = 1*10**(-itr)

    for dim_itr in range(2):

        print("\n\n##########################################################################")
        print(">> Starting finite difference level ",itr, " dimension ",dim_itr+1)
        print("##########################################################################")

        # Create case folder
        case_dir = "tmp_calculation_"+str(itr)+"_"+str(dim_itr+1)
        if os.path.exists(case_dir):
            os.system("rm -rf "+case_dir+"/")
        os.system("mkdir "+case_dir)

        # Copy all necessary files to case folder
        shutil.copy("ubend_2d.mdpa", case_dir+"/")
        shutil.copy("OPT_model.mdpa", case_dir+"/")
        shutil.copy("parameters_optimization.json", case_dir+"/")
        shutil.copy("restart_flow.dat", case_dir+"/")
        shutil.copy("ubend_coarse.cfg", case_dir+"/")
        shutil.copy("ubend_2d.su2", case_dir+"/")
        shutil.copy("su2_mesh_data.json", case_dir+"/")

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

        OutputOPTAsGid(opt_part, "opt_part")

        # Update mesh
        kso.MeshControllerUtilities(opt_part).UpdateMeshAccordingInputVariable(kso.SHAPE_UPDATE)
        kso.MeshControllerUtilities(opt_part).SetReferenceMeshToMesh()

        # Run analysis
        analyzer.InitializeBeforeOptimizationLoop()
        cfd_response_value = analyzer.AnalyzeDesignAndReportToCommunicator(opt_part, 1)

        # Evaluate fd gradients
        if dim_itr == 0:
            cfd_grad_x = (cfd_response_value - cfd_response_value_reference) / current_delta
        elif dim_itr == 1:
            cfd_grad_y = (cfd_response_value - cfd_response_value_reference) / current_delta

        print("\n> current cfd response value =",cfd_response_value)
        print("\n")

        # Change back to top dir
        os.chdir(top_dir)

    # Write results
    with open("fd_results_"+str(response_function)+"_move_node_"+str(move_node_id)+".txt",'a') as results_file:
        results_file.write(str(current_delta)+","+str(cfd_grad_x)+","+str(cfd_grad_y)+","+str(time.ctime())+"\n")

print("\n\n##########################################################################")
print("\n> Finished verifiation process in" ,round( time.time()-start_time, 3 ), " s.")
print("##########################################################################")
