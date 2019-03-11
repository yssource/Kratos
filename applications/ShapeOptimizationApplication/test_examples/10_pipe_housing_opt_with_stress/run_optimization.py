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

# Read parameters
with open("parameters_optimization.json",'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())

optimization_model = km.Model()

# Definition of external analyzer
class CustomAnalyzer(AnalyzerBaseClass):
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __init__(self):
        self.cfd_part_name = "CFD_part"
        self.cfd_part = None
        self.cfd_interface_part_name = "design_surface"
        self.cfd_damping_part_name = "edge_nodes"
        self.cfd_interface_part = None
        self.cfd_interface_damping_utils = None

        self.optimization_mdpa_name = "CSM_domain"
        self.csm_interface_part_name = "design_surface_wet_interface"
        self.optimization_part = None

        self.material_filename = "parameters_material.json"
        self.results_folder = "Optimization_Results"
        self.reference_pressure = 74571.25

        self.cfd_interface_gid_output = None
        self.cfd_interface_output_parameters = km.Parameters("""
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
                "nodal_results"       : ["TRACTION", "NORMAL", "DF1DX", "MESH_CHANGE", "MESH_DISPLACEMENT"],
                "gauss_point_results" : []
            },
            "point_data_configuration"  : []
        }""")
        self.adjoint_output_parameters = km.Parameters("""
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
                "nodal_results"       : ["ADJOINT_DISPLACEMENT","SHAPE_SENSITIVITY"],
                "gauss_point_results" : []
            },
            "point_data_configuration"  : []
        }""")
        self.traction_mapper_settings = km.Parameters("""
        {
            "mapper_type" : "nearest_element",
            "echo_level"  : 0
        }""")

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        # Initialize elements to avoid seg-fault when not using any internal analyzer but writing the whole optimization model part with gid (which contains elements)
        self.optimization_part = optimization_model.GetModelPart(self.optimization_mdpa_name)
        material_settings = km.Parameters("""{"Parameters": {"materials_filename": ""}} """)
        material_settings["Parameters"]["materials_filename"].SetString(self.material_filename)
        km.ReadMaterialsUtility(material_settings, optimization_model)
        for elem in self.optimization_part.Elements:
            elem.Initialize()

        # Initialize cfd data
        cfd_model = km.Model()
        self.cfd_part = cfd_model.CreateModelPart(self.cfd_part_name)
        self.cfd_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)
        self.cfd_part.AddNodalSolutionStepVariable(kso.TRACTION)
        self.cfd_part.AddNodalSolutionStepVariable(kso.DF1DX)
        self.cfd_part.AddNodalSolutionStepVariable(km.MESH_DISPLACEMENT)
        self.cfd_part.AddNodalSolutionStepVariable(kso.MESH_CHANGE)
        self.cfd_part.AddNodalSolutionStepVariable(km.NORMAL)
        self.cfd_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)

        model_part_io = km.ModelPartIO(self.cfd_part_name)
        model_part_io.ReadModelPart(self.cfd_part)

        self.cfd_interface_part = self.cfd_part.GetSubModelPart(self.cfd_interface_part_name)

        # Initialize damping of mesh displacement on CFD interface
        damping_regions = {}
        damping_regions[self.cfd_damping_part_name] = self.cfd_part.GetSubModelPart(self.cfd_damping_part_name)
        modified_settings_for_damping = parameters["optimization_settings"]["model_settings"]["damping"].Clone()
        modified_settings_for_damping["damping_regions"][0]["sub_model_part_name"].SetString(self.cfd_damping_part_name)
        self.cfd_interface_damping_utils = kso.DampingUtilities(self.cfd_part, damping_regions, modified_settings_for_damping)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        # get mesh motion
        mesh_change = {}
        for node in self.optimization_part.Nodes:
            mesh_change[node.Id] = node.GetSolutionStepValue(kso.MESH_CHANGE)

        # Update fluid mesh (structure is controlled by the optimization algorithm)
        vm_mapper = kso.MapperVertexMorphingMatrixFree(current_design, self.cfd_interface_part, parameters["optimization_settings"]["design_variables"]["filter"])
        vm_mapper.Map(kso.CONTROL_POINT_UPDATE, km.MESH_DISPLACEMENT)

        if parameters["optimization_settings"]["model_settings"]["damping"]["apply_damping"].GetBool():
            self.cfd_interface_damping_utils.DampNodalVariable(km.MESH_DISPLACEMENT)

        kso.MeshControllerUtilities(self.cfd_interface_part).UpdateMeshAccordingInputVariable(km.MESH_DISPLACEMENT)
        kso.MeshControllerUtilities(self.cfd_interface_part).SetReferenceMeshToMesh()
        kso.MeshControllerUtilities(self.cfd_interface_part).LogMeshChangeAccordingInputVariable(km.MESH_DISPLACEMENT)
        kso.GeometryUtilities(self.cfd_interface_part).ComputeUnitSurfaceNormals()

        # Run fluid

        # Read and dimensionalize fluid results
        nodal_force_densities = self.__ReadNodalValuesFromCSVFile("surface_flow_DSN001.csv",0,[9,10,11],1,',')
        for node in self.cfd_interface_part.Nodes:
            kratos_node_id = node.Id
            su2_node_id = kratos_node_id
            if kratos_node_id not in nodal_force_densities.keys():
                su2_node_id = 0
                print("WARNING!!!!! Node-Id found in SU2 Ids which do not have a coutnerpart in Kratos. Assuming this corresponds to the SU2 node with a 0 node-Id.")
            traction = nodal_force_densities[su2_node_id]
            dimensional_traction = [value*self.reference_pressure for value in traction]
            # # To consider another pressure level
            # new_pressure_level = 3e7 # 1e5 Pa --> 1 bar
            # dimensional_traction = np.array(dimensional_traction) + new_pressure_level * np.array(node.GetSolutionStepValue(kso.NORMALIZED_SURFACE_NORMAL))
            # dimensional_traction = dimensional_traction.tolist()
            node.SetSolutionStepValue(kso.TRACTION, dimensional_traction)

        # Output results on CFD interface
        self.__OutputCFDInterface(optimization_iteration)

        # Solve CSM
        displacements, value, gradients = self.__RunCSM(partition_itr=1, opt_itr=optimization_iteration, mesh_displacement=mesh_change, displacements={}, cfd_interface_part=self.cfd_interface_part)

        if communicator.isRequestingValueOf("mean_stress_partition_1"):
            communicator.reportValue("mean_stress_partition_1", value)

        if communicator.isRequestingGradientOf("mean_stress_partition_1"):
            communicator.reportGradient("mean_stress_partition_1", gradients)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __RunCSM(self, partition_itr, opt_itr, mesh_displacement, displacements, cfd_interface_part):
        print("\n> Starting __RunCSM in analyzer...")
        start_time = time.time()

        analysis_model = km.Model()

        response_settings = parameters["optimization_settings"]["objectives"][0]["kratos_response_settings"]
        response_settings["critical_part_name"].SetString("stress_partition_"+str(partition_itr))
        response_settings["primal_settings"].SetString("parameters_analysis_partition_"+str(partition_itr)+".json")

        csm = structural_response_function_factory.CreateResponseFunction(response_settings["response_type"].GetString(), response_settings, analysis_model)

        csm_primal_part = csm.primal_model_part
        csm_primal_part.AddNodalSolutionStepVariable(kso.MESH_CHANGE)
        csm_primal_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)
        csm_primal_part.AddNodalSolutionStepVariable(km.NORMAL)
        csm_primal_part.AddNodalSolutionStepVariable(kso.TRACTION)
        csm_adjoint_part = csm.adjoint_model_part

        # Initialize and set time to current optimization iteration
        csm.primal_analysis.project_parameters["problem_data"]["start_time"].SetDouble(opt_itr-1)
        csm.primal_analysis.project_parameters["problem_data"]["end_time"].SetDouble(opt_itr)

        csm.adjoint_analysis.project_parameters["problem_data"]["start_time"].SetDouble(opt_itr-1)
        csm.adjoint_analysis.project_parameters["problem_data"]["end_time"].SetDouble(opt_itr)

        csm.Initialize()

        csm_primal_part.ProcessInfo.SetValue(km.STEP, opt_itr-1)
        csm_adjoint_part.ProcessInfo.SetValue(km.STEP, opt_itr-1)

        # apply mesh motion
        if mesh_displacement != {}:
            for node in csm_primal_part.Nodes:
                mesh_change = mesh_displacement[node.Id]
                node.SetSolutionStepValue(kso.MESH_CHANGE, mesh_change)
                node.X += mesh_change[0]
                node.Y += mesh_change[1]
                node.Z += mesh_change[2]
                node.X0 += mesh_change[0]
                node.Y0 += mesh_change[1]
                node.Z0 += mesh_change[2]

        # compute primal field if necessary
        if displacements == {}:

            # Map forces
            csm_interface_part = csm_primal_part.GetSubModelPart(self.csm_interface_part_name)
            kso.GeometryUtilities(csm_interface_part).ComputeUnitSurfaceNormals()

            cfd_to_csm_mapper = kma.MapperFactory.CreateMapper(cfd_interface_part, csm_interface_part, self.traction_mapper_settings.Clone())
            cfd_to_csm_mapper.Map(kso.TRACTION, kso.TRACTION)

            # Apply forces
            for node_i in csm_interface_part.Nodes:
                normal = node_i.GetSolutionStepValue(km.NORMAL)
                area = math.sqrt( normal[0]**2+normal[1]**2+normal[2]**2 )
                traction = node_i.GetSolutionStepValue(kso.TRACTION)
                force = [value*area for value in traction]
                node_i.SetSolutionStepValue(kso.TRACTION, traction)
                node_i.SetSolutionStepValue(kcsm.POINT_LOAD, force)

            # Compute primals
            csm.InitializeSolutionStep()

            # Store primal result
            displacements = {node.Id: node.GetSolutionStepValue(km.DISPLACEMENT) for node in csm_primal_part.Nodes}
        else:
            csm_analys = csm.primal_analysis
            csm_analys._GetSolver().AdvanceInTime(current_time=opt_itr-1)
            csm_analys.InitializeSolutionStep()

            for node in csm_primal_part.Nodes:
                node.SetSolutionStepValue(km.DISPLACEMENT, displacements[node.Id])

            csm_analys.FinalizeSolutionStep()
            csm_analys.OutputSolutionStep()

        csm.CalculateValue()
        csm.CalculateGradient()
        csm.FinalizeSolutionStep()
        csm.Finalize()

        # Some postprocessing
        primal_results_filename = "primal_results_partition_"+str(partition_itr)+".post.bin"
        adjoint_results_filename = "adjoint_results_partition_"+str(partition_itr)+".post.bin"
        self.__OutputMdpaAsGid(csm.adjoint_model_part, adjoint_results_filename.replace(".post.bin",""), self.adjoint_output_parameters)
        self.__CopyBinFileToResultsFolder(primal_results_filename, opt_itr)
        self.__CopyBinFileToResultsFolder(adjoint_results_filename, opt_itr)

        print("> Finished __RunCSM in" ,round( time.time()-start_time, 3 ), " s.")

        return displacements, csm.GetValue(), csm.GetShapeGradient()

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __OutputCFDInterface( self, optimization_iteration ):
        self.cfd_interface_part.ProcessInfo[km.TIME] = optimization_iteration
        if self.cfd_interface_gid_output is None:
            self.cfd_interface_gid_output = GiDOutputProcess( self.cfd_interface_part, os.path.join(self.results_folder,"CFD_interface"), self.cfd_interface_output_parameters)
            self.cfd_interface_gid_output.ExecuteInitialize()
            self.cfd_interface_gid_output.ExecuteBeforeSolutionLoop()
        self.cfd_interface_gid_output.ExecuteInitializeSolutionStep()
        self.cfd_interface_gid_output.PrintOutput()
        self.cfd_interface_gid_output.ExecuteFinalizeSolutionStep()

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __CopyBinFileToResultsFolder(self, filename, optimization_iteration):
        old_name = filename
        new_name = old_name.replace(".post.bin","_DESIGN_"+str(optimization_iteration)+".post.bin")
        os.rename(old_name,new_name)
        shutil.move(new_name, self.results_folder)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def __OutputMdpaAsGid( output_mdpa, output_filename, parameters ):
        gid_output_original = GiDOutputProcess( output_mdpa, output_filename, parameters )
        gid_output_original.ExecuteInitialize()
        gid_output_original.ExecuteBeforeSolutionLoop()
        gid_output_original.ExecuteInitializeSolutionStep()
        gid_output_original.PrintOutput()
        gid_output_original.ExecuteFinalizeSolutionStep()
        gid_output_original.ExecuteFinalize()

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def __ReadNodalValuesFromCSVFile(filename, collum_with_node_id, collums_to_read, number_of_header_rows, delimiter_sign):
        values = {}

        if len(collums_to_read) == 1:
            with open(filename, newline='') as csvfile:
                csv_reader = csv.reader(csvfile, delimiter=delimiter_sign, quotechar='|')
                for itr in range(number_of_header_rows):
                    csvfile.readline()
                for row in csv_reader:
                    node_id = int(row[collum_with_node_id])
                    values[node_id] = float(row[collums_to_read[0]])
        else:
            with open(filename, newline='') as csvfile:
                csv_reader = csv.reader(csvfile, delimiter=delimiter_sign, quotechar='|')
                for itr in range(number_of_header_rows):
                    csvfile.readline()
                for row in csv_reader:
                    node_id = int(row[collum_with_node_id])
                    entries = []
                    for coll_id in collums_to_read:
                        entries.append(float(row[coll_id]))
                    values[node_id] = entries

        return values


    # -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model, CustomAnalyzer())
optimizer.Optimize()