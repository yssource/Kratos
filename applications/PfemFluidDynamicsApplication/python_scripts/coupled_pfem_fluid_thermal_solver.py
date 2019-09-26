from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledPfemFluidThermalSolver(main_model_part, custom_settings)

class CoupledPfemFluidThermalSolver(PythonSolver):

    def __init__(self, model, custom_settings):

        super(CoupledPfemFluidThermalSolver, self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "coupled_pfem_fluid_thermal_solver",
            "model_part_name": "PfemFluidModelPart",
            "echo_level"                         : 1,
            "fluid_solver_settings":{
                "physics_type"   : "fluid",
                "domain_size": 2,
                "time_stepping"               : {
                    "automatic_time_step" : false,
                    "time_step"           : 0.001
                },
                "model_import_settings":{
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                },
                "buffer_size": 3,
                "echo_level": 1,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": true,
                "move_mesh_flag": true,
                "dofs"                : [],
                "stabilization_factor": 1.0,
                "line_search": false,
                "compute_contact_forces": false,
                "block_builder": false,
                "component_wise": false,
                "predictor_corrector": true,
                "time_order": 2,
                "maximum_velocity_iterations": 1,
                "maximum_pressure_iterations": 7,
                "velocity_tolerance": 1e-5,
                "pressure_tolerance": 1e-5,
                "pressure_linear_solver_settings":  {
                    "solver_type"                    : "amgcl",
                    "max_iteration"                  : 5000,
                    "tolerance"                      : 1e-9,
                    "provide_coordinates"            : false,
                    "scaling"                        : false,
                    "smoother_type"                  : "damped_jacobi",
                    "krylov_type"                    : "cg",
                    "coarsening_type"                : "aggregation",
                    "verbosity"                      : 0
                },
                "velocity_linear_solver_settings": {
                    "solver_type"                    : "bicgstab",
                    "max_iteration"                  : 5000,
                    "tolerance"                      : 1e-9,
                    "preconditioner_type"            : "none",
                    "scaling"                        : false
                },
                "solving_strategy_settings":{
                   "time_step_prediction_level": 0,
                   "max_delta_time": 1.0e-5,
                   "fraction_delta_time": 0.9,
                   "rayleigh_damping": false,
                   "rayleigh_alpha": 0.0,
                   "rayleigh_beta" : 0.0
                },
                "bodies_list": [],
                "problem_domain_sub_model_part_list": [],
                "processes_sub_model_part_list": [],
                "constraints_process_list": [],
                "loads_process_list"       : [],
                "output_process_list"      : [],
                "output_configuration"     : {},
                "problem_process_list"     : [],
                "processes"                : {},
                "output_processes"         : {},
                "check_process_list": []
            },
            "thermal_solver_settings": {
                "solver_type": "Transient",
                "analysis_type": "linear",
                "computing_model_part_name": "fluid_computing_domain",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                }
            }
        }
        """)
            
        #self.settings = custom_settings
        ## Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Get domain size
        self.domain_size = self.settings["fluid_solver_settings"]["domain_size"].GetInt()
        
        from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_fluid_solver
        self.fluid_solver = pfem_fluid_solver.CreateSolver(self.model,self.settings["fluid_solver_settings"]) 

        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")

    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both fluid and thermal solvers.
        self.fluid_solver.AddVariables()
        self.thermal_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)
        print("::[Coupled Pfem Fluid Thermal Solver]:: Variables MERGED")
        ## 2190920: a warning is raised when the thermal model is cloned form the fluid one because the pfem extra variables are
        ##          added after (by the pfem_analysis and not by pfem_fluid_solver)

    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa
        self.fluid_solver.ImportModelPart() # import model fluid model part and call pfem_check_and_prepare_model_process_fluid
        # self.fluid_solver._ImportModelPart(self.fluid_solver.main_model_part,self.settings["fluid_solver_settings"]["model_import_settings"])

        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
            #modeler.GenerateModelPart(self.fluid_solver.GetComputingModelPart(),
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff2D",
                                      "ThermalFace2D2N")
        else:
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff3D",
                                      "ThermalFace3D3N")
        print("::[Coupled Pfem Fluid Thermal Solver]:: Thermal_model_part CLONED")
        
        # self.fluid_solver.ImportModelPart()

        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

        # self.settings["thermal_solver_settings"]["computing_model_part_name"].SetString("fluid_computing_domain")
        print(1)

    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        self.thermal_solver.PrepareModelPart()

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()

    def InitializeStrategy(self):
        self.fluid_solver.InitializeStrategy()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        self.thermal_solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def Solve(self):
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()
