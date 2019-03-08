# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
import KratosMultiphysics.StructuralMechanicsApplication as KCSM
import structural_response_function_factory

# Additional imports
from analyzer_base import AnalyzerBaseClass
from gid_output_process import GiDOutputProcess
import time

# Read parameters
with open("parameters_optimization.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

optimization_model = KM.Model()

# Definition of external analyzer
class CustomAnalyzer(AnalyzerBaseClass):
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __init__(self):
        self.optimization_mdpa_name = "CSM_domain"
        self.optimization_part = None
        self.wet_interface_part_name_csm = "design_surface_wet_interface"
        self.material_filename = "parameters_material.json"
        self.adjoint_output_parameters = KM.Parameters("""
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

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        # Initialize elements to avoid seg-fault when not using any internal analyzer but writing the whole optimization model part with gid (which contains elements)
        self.optimization_part = optimization_model.GetModelPart(self.optimization_mdpa_name)
        material_settings = KM.Parameters("""{"Parameters": {"materials_filename": ""}} """)
        material_settings["Parameters"]["materials_filename"].SetString(self.material_filename)
        KM.ReadMaterialsUtility(material_settings, optimization_model)
        for elem in self.optimization_part.Elements:
            elem.Initialize()

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        # Determine nodal forces
        interface_part_opt = self.optimization_part.GetSubModelPart(self.wet_interface_part_name_csm)
        KSO.GeometryUtilities(interface_part_opt).ComputeUnitSurfaceNormals()

        interface_forces = {}
        pressure = 1e5
        for node in interface_part_opt.Nodes:
            normal = node.GetSolutionStepValue(KM.NORMAL)
            interface_forces[node.Id] = [-value*pressure for value in normal]

        # Solve stress different stress partitions
        displacements, gradients = self.__RunCSM(itr=1, mesh_displacement={}, displacements={}, interface_forces=interface_forces)
        _, gradients = self.__RunCSM(itr=2, mesh_displacement={}, displacements=displacements, interface_forces={})
        _, gradients = self.__RunCSM(itr=3, mesh_displacement={}, displacements=displacements, interface_forces={})
        _, gradients = self.__RunCSM(itr=4, mesh_displacement={}, displacements=displacements, interface_forces={})

        # # Report values
        # if communicator.isRequestingValueOf("lifetime"):
        #     communicator.reportValue("lifetime", value)

        # if communicator.isRequestingGradientOf("lifetime"):
        #     response.CalculateGradient()
        #     communicator.reportGradient("lifetime", response.GetShapeGradient())


    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __RunCSM(self, itr, mesh_displacement, displacements, interface_forces):
        print("\n> Starting __RunCSM in analyzer...")
        start_time = time.time()

        analysis_model = KM.Model()

        response_settings = parameters["optimization_settings"]["objectives"][0]["kratos_response_settings"]
        response_settings["critical_part_name"].SetString("stress_partition_"+str(itr))
        response_settings["primal_settings"].SetString("parameters_analysis_partition_"+str(itr)+".json")

        csm = structural_response_function_factory.CreateResponseFunction(response_settings["response_type"].GetString(), response_settings, analysis_model)

        csm_primal_part = csm.primal_model_part
        csm_primal_part.AddNodalSolutionStepVariable(KSO.MESH_CHANGE)
        csm_primal_part.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)
        csm_primal_part.AddNodalSolutionStepVariable(KM.NORMAL)

        csm.Initialize()

        # apply mesh motion
        if mesh_displacement != {}:
            pass

        # compute primal field if necessary
        if displacements == {}:

            # Apply forces
            interface_part = csm_primal_part.GetSubModelPart(self.wet_interface_part_name_csm)
            for node_i in interface_part.Nodes:
                node_i.SetSolutionStepValue(KCSM.POINT_LOAD, interface_forces[node_i.Id])

            # Compute primals
            csm.InitializeSolutionStep()

            # Store primal result
            displacements = {node.Id: node.GetSolutionStepValue(KM.DISPLACEMENT) for node in csm_primal_part.Nodes}
        else:
            csm_analys = csm.primal_analysis
            csm_analys._GetSolver().AdvanceInTime(current_time=0)
            csm_analys.InitializeSolutionStep()

            for node in csm_primal_part.Nodes:
                node.SetSolutionStepValue(KM.DISPLACEMENT, displacements[node.Id])

            csm_analys.FinalizeSolutionStep()
            csm_analys.OutputSolutionStep()

        csm.CalculateValue()
        csm.CalculateGradient()
        csm.FinalizeSolutionStep()
        csm.Finalize()

        self.__OutputMdpaAsGid(csm.adjoint_model_part, "adjoint_results_csm_"+str(itr), self.adjoint_output_parameters)

        print("> Finished __RunCSM in" ,round( time.time()-start_time, 3 ), " s.")

        return displacements, csm.GetShapeGradient()

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

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model, CustomAnalyzer())
optimizer.Optimize()