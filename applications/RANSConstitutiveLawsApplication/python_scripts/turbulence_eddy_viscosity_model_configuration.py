import KratosMultiphysics as KM
import KM.RANSConstitutiveLawsApplication as KratosRANS
from KM.kratos_utilities import IsApplicationAvailable

if IsApplicationAvailable("FluidDynamicsApplication"):
    import KM.FluidDynamicsApplication as KratosCFD
    from turbulence_model_configuration import TurbulenceModelConfiguration
else:
    raise Exception("RANSConstitutiveLawsApplication requires FluidDynamicsApplication which is not found. Please install/compile it and try again.")

class RANSSolverConfigurations(object):
    def __init__(self):
        self.elements = []
        self.conditions = []
        self.linear_solvers = []
        self.time_schemes = []
        self.model_parts = []
        self.conv_criterias = []
        self.builder_and_solvers = []
    
    def AddSolver(self, element_name):
        pass

class TurbulenceEddyViscosityModelConfiguration(TurbulenceModelConfiguration):
    def __init__(self, model, settings):
        default_settings = KM.Parameters(r'''{
            "model_type"            : "",
            "fluid_model_part"      : "PLEASE_SPECIFY_FLUID_MODEL_PART",
            "inlet_conditions"      : ["PLEASE_SPECIFY_INLET_CONDITIONS"],
            "outlet_conditions"     : ["PLEASE_SPECIFY_OUTLET_CONDITIONS"],
            "wall_conditions"       : ["PLEASE_SPECIFY_WALL_CONDITIONS"],
            "distance_calculation"  : { 
                "max_iterations"         : 5,
                "linear_solver_settings" : {}
            },
            "y_plus_calculation"    : {
                "model_type"        : ""
                "model_settings"    : {}
            }
            "mesh_moving"       : false,
            "echo_level"        : 0,
            "model_settings"    : {},
            "coupling_settings" : {
                "relative_tolerance"    : 1e-3,
                "absolute_tolerance"    : 1e-5,
                "max_iterations"        : 200,
                "echo_level"            : 0                    
            }
        }''')

        self.settings = settings.ValidateAndAssignDefaults(default_settings)

        super(TurbulenceEddyViscosityModelConfiguration, self).__init__(model, settings)

        self.fluid_model_part = self.model.GetSubModelPart(
            self.settings["fluid_model_part"].GetString())

        self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.distance_calculation_process = None
        self.y_plus_model_process = None
        self.strategies_list = []

    def Initialize(self):
        turbulence_model_parts_list = self.CreateTurbulenceModelParts(self.fluid_model_part)
        
        self.__InitializeNodeFlags()        
        self.__CalculateWallDistances()

        self.GetTurbulenceModelProcess().ExecuteInitialize()
   
    def InitializeSolutionStep(self):
        if self.mesh_moving:
            self.__InitializeNodeFlags()        
            self.__CalculateWallDistances()

        self.GetTurbulenceModelProcess().ExecuteInitializeSolutionStep()
    
    def SolveSolutionStep(self):
        for strategy in self.strategies_list:
            strategy.InitializeSolutionStep()
            strategy.Predict()

        is_converged = False
        iteration = 1

        variable_utils = KM.VariableUtils()

        while (not is_converged and iteration <= self.settings["coupling_settings"]["max_iterations"].GetInt()):
            variable_utils.CopyScalarVar(KM.OLD_TURBULENT_VISCOSITY, KM.TURBULENT_VISCOSITY, self.fluid_model_part.Nodes)                    
            for strategy in self.strategies_list:
                strategy.SolveSolutionStep()

            self.UpdateTurbulentViscosity()

            iteration += 1
                    
        for strategy in self.strategies_list:
            strategy.FinalizeSolutionStep()
        
        self.UpdateFluidViscosity()

    def FinalizeSolutionStep(self):
        self.GetTurbulenceModelProcess().ExecuteFinalizeSolutionStep()
        self.UpdateFluidViscosity()

    def AddVariables(self):
        self.fluid_model_part.AddNodalSolutionStepVariable(RANS_Y_PLUS)        
        self.fluid_model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(KM.FLAG_VARIABLE)
        self.fluid_model_part.AddNodalSolutionStepVariable(KM.KINEMATIC_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(KM.TURBULENT_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(KM.OLD_TURBULENT_VISCOSITY)

    def AddDofs(self):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::AddDofs"
        )

    def GetTurbulenceModelProcess(self):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::GetTurbulenceModelProcess"
        )
    
    def CalculateYPlus(self):
        if self.y_plus_model_process is None:
            import rans_y_plus_models
            self.y_plus_model_process = rans_y_plus_models.Factory(self.fluid_model_part, self.settings["y_plus_calculation"])
            KM.Logger(self.__name__, "Initialized " + self.y_plus_model_process)
        
        self.y_plus_model_process.Execute()
        KM.Logger(self.__name__, "Calculated y plus using " + self.y_plus_model_process)

    def UpdateTurbulentViscosity(self):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::UpdateTurbulentViscosity"
        )        

    def UpdateFluidViscosity(self):
        rans_variable_utils = KratosRANS.RansVariableUtils()
        rans_variable_utils.AddToHistoricalNodeScalarVariable(KM.VISCOSITY, KM.KINEMATIC_VISCOSITY, KM.TURBULENT_VISCOSITY, self.fluid_model_part)
    
    def __InitializeNodeFlags(self):
        self.__InitializeNodeFlags(self.settings["inlet_conditions"].GetStringArray(), KM.INLET, True)
        self.__InitializeNodeFlags(self.settings["outlet_conditions"].GetStringArray(), KM.OUTLET, True)
        self.__InitializeNodeFlags(self.settings["wall_conditions"].GetStringArray(), KM.STRUCTURE, True)

        variable_utils = KM.VariableUtils()
        variable_utils.SetScalarVar(KM.DISTANCE, 1.0, self.fluid_model_part.Nodes)
        variable_utils.SetScalarVar(KM.DISTANCE, 0.0, self.fluid_model_part.Nodes, KM.STRUCTURE)
        
        variable_utils.CopyScalarVar(KM.VISCOSITY, KM.KINEMATIC_VISCOSITY, self.fluid_model_part.Nodes)        
    
    def __InitializeNodeFlags(self, conditions_list, flag, value):
        variable_utils = KM.VariableUtils()
        for condition_name in conditions_list:
            if not self.fluid_model_part.HasSubModelPart(condition_name):
                raise Exception(condition_name + " not found in " + self.fluid_model_part)
            variable_utils.SetFlag(flag, value, self.GetSubModelPart(condition_name).Nodes)

    def __CalculateWallDistances(self):
        if (self.distance_calculation_process is None):
            import KM.python_linear_solver_factory as linear_solver_factory
            self.distance_calculation_linear_solver = linear_solver_factory.ConstructSolver(self.settings["distance_calculation"]["linear_solver_settings"])
            max_iterations = self.settings["distance_calculation"]["max_iterations"].GetInt()
            if (self.domain_size == 2):
                self.distance_calculation_process = KM.VariationalDistanceCalculationProcess2D(self.fluid_model_part, self.distance_calculation_linear_solver, max_iterations)
            elif (self.domain_size == 3):
                self.distance_calculation_process = KM.VariationalDistanceCalculationProcess3D(self.fluid_model_part, self.distance_calculation_linear_solver, max_iterations)
            else:
                raise Exception("Unsupported domain size")
            
            KM.Logger.PrintInfo(self.__name__, "Variational distance calculation process initialized.")
        
        self.distance_calculation_process.Execute()
        KM.Logger(self.__name__, "Wall distances calculated.")
