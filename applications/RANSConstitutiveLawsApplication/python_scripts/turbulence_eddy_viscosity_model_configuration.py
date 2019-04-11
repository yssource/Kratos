import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSConstitutiveLawsApplication as KratosRANS
import math

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
    from turbulence_model_configuration import TurbulenceModelConfiguration
else:
    raise Exception(
        "RANSConstitutiveLawsApplication requires FluidDynamicsApplication which is not found. Please install/compile it and try again."
    )

class TurbulenceEddyViscosityModelConfiguration(TurbulenceModelConfiguration):
    def __init__(self, model, settings):
        default_settings = Kratos.Parameters(r'''{
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
                "model_type"        : "",
                "model_settings"    : {}
            },
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

        super(TurbulenceEddyViscosityModelConfiguration, self).__init__(
            model, settings)

        self.fluid_model_part = self.model.GetModelPart(
            self.settings["fluid_model_part"].GetString())

        self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.distance_calculation_process = None
        self.y_plus_model_process = None
        self.strategies_list = []

        self.relative_tolerance = self.settings["coupling_settings"][
            "relative_tolerance"].GetDouble()
        self.absolute_tolerance = self.settings["coupling_settings"][
            "absolute_tolerance"].GetDouble()

        self.is_computing_solution = False

    def Initialize(self):
        for model_part in self.turbulence_model_parts_list:
            self.AddDofs(model_part)

        self.PrepareSolvingStrategy()

        for strategy in self.strategies_list:
            strategy.Initialize()

        self.__InitializeModelPart()

        variable_utils = Kratos.VariableUtils()
        variable_utils.CopyScalarVar(Kratos.VISCOSITY, Kratos.KINEMATIC_VISCOSITY,  self.fluid_model_part.Nodes)

    def InitializeSolutionStep(self):
        if self.mesh_moving:
            self.__InitializeModelPart()

    def Execute(self):
        print("help trampoline")
        if not self.is_computing_solution:
            return

        self.UpdateBoundaryConditions()

        for strategy in self.strategies_list:
            strategy.InitializeSolutionStep()
            strategy.Predict()

        is_converged = False
        iteration = 1

        variable_utils = Kratos.VariableUtils()
        while (not is_converged and iteration <=
               self.settings["coupling_settings"]["max_iterations"].GetInt()):
            variable_utils.CopyScalarVar(KratosRANS.OLD_TURBULENT_VISCOSITY,
                                         KratosCFD.TURBULENT_VISCOSITY,
                                         self.fluid_model_part.Nodes)

            # solve for turbulent scalar variables
            for strategy in self.strategies_list:
                strategy.SolveSolutionStep()

            # calculating the new turbulent viscosity
            self.UpdateTurbulentViscosity()

            # checking for convergence of turbulent viscosity
            increase_norm = KratosRANS.RansVariableUtils(
            ).GetScalarVariableIncreaseNormSquare(
                self.fluid_model_part.Nodes,
                KratosRANS.OLD_TURBULENT_VISCOSITY, Kratos.TURBULENT_VISCOSITY)

            solution_norm = KratosRANS.RansVariableUtils(
            ).GetScalarVariableSolutionNormSquare(self.fluid_model_part.Nodes,
                                                  Kratos.TURBULENT_VISCOSITY)
            if (solution_norm == 0.0):
                solution_norm = 1.0

            relative_tolerance = increase_norm / solution_norm
            absolute_tolerance = math.sqrt(
                increase_norm) / self.fluid_model_part.NumberOfNodes()

            is_converged = (relative_tolerance < self.relative_tolerance
                            or absolute_tolerance < self.absolute_tolerance)
            iteration += 1

        for strategy in self.strategies_list:
            strategy.FinalizeSolutionStep()

        self.UpdateFluidViscosity()

    def AddVariables(self, model_part):
        # adding variables from navier_stokes_vms_monolithic_solver
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.IS_STRUCTURE)
        model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        model_part.AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        model_part.AddNodalSolutionStepVariable(Kratos.NODAL_H)
        model_part.AddNodalSolutionStepVariable(Kratos.ADVPROJ)
        model_part.AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        model_part.AddNodalSolutionStepVariable(Kratos.Y_WALL)
        model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        # adding variables required by rans
        model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)
        model_part.AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        model_part.AddNodalSolutionStepVariable(KratosRANS.OLD_TURBULENT_VISCOSITY)

    def AddDofs(self, model_part):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::AddDofs"
        )

    def UpdateTurbulentViscosity(self):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::UpdateTurbulentViscosity"
        )

    def UpdateBoundaryConditions(self):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::UpdateBoundaryConditions"
        )

    def InitializeBoundaryNodes(self):
        raise Exception(
            "Calling base class TurbulenceEddyViscosityModelConfiguration::InitializeBoundaryNodes"
        )

    def CalculateYPlus(self):
        if self.y_plus_model_process is None:
            import rans_y_plus_models
            self.y_plus_model_process = rans_y_plus_models.Factory(
                self.fluid_model_part, self.settings["y_plus_calculation"])
            Kratos.Logger(self.__name__,
                      "Initialized " + self.y_plus_model_process)

        self.y_plus_model_process.Execute()
        Kratos.Logger(self.__name__,
                  "Calculated y plus using " + self.y_plus_model_process)

    def UpdateFluidViscosity(self):
        rans_variable_utils = KratosRANS.RansVariableUtils()
        rans_variable_utils.AddToHistoricalNodeScalarVariable(
            Kratos.VISCOSITY, Kratos.KINEMATIC_VISCOSITY, Kratos.TURBULENT_VISCOSITY,
            self.fluid_model_part)

    def CreateStrategy(self, solver_settings, scheme_settings, model_part,
                       scalar_variable, scalar_variable_rate,
                       relaxed_scalar_variable_rate):
        import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            solver_settings["linear_solver_settings"])
        convergence_criteria = KratosRANS.GenericScalarConvergenceCriteria(
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())
        builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(
            linear_solver)
        time_scheme = KratosRANS.GenericResidualBasedBossakVelocityDynamicScalarScheme(
            scheme_settings["alpha_bossak"].GetDouble(), scalar_variable,
            scalar_variable_rate, relaxed_scalar_variable_rate)

        strategy = Kratos.ResidualBasedNewtonRaphsonStrategy(
            model_part, time_scheme, linear_solver, convergence_criteria,
            builder_and_solver, solver_settings["max_iterations"].GetInt(),
            False, False, False)

        return strategy, linear_solver, convergence_criteria, builder_and_solver, time_scheme

    def __InitializeModelPart(self):
        self.__InitializeNodeFlags()
        self.InitializeBoundaryNodes()

        variable_utils = Kratos.VariableUtils()
        variable_utils.SetScalarVar(Kratos.DISTANCE, 1.0, self.fluid_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.DISTANCE, 0.0, self.fluid_model_part.Nodes, Kratos.STRUCTURE, True)

        self.__CalculateWallDistances()

    def __InitializeNodeFlags(self):
        self.__InitializeNodeFlagsForConditions(
            self.settings["inlet_conditions"].GetStringArray(), Kratos.INLET, True)
        self.__InitializeNodeFlagsForConditions(self.settings["outlet_conditions"].GetStringArray(), Kratos.OUTLET, True)
        self.__InitializeNodeFlagsForConditions(self.settings["wall_conditions"].GetStringArray(), Kratos.STRUCTURE, True)
        self.__InitializeNodeFlagsForConditions(self.settings["wall_conditions"].GetStringArray(), Kratos.INLET, False)
        self.__InitializeNodeFlagsForConditions(self.settings["wall_conditions"].GetStringArray(), Kratos.INLET, False)


    def __InitializeNodeFlagsForConditions(self, conditions_list, flag, value):
        variable_utils = Kratos.VariableUtils()
        for condition_name in conditions_list:
            if not self.fluid_model_part.HasSubModelPart(condition_name):
                raise Exception(condition_name + " not found in " +
                                self.fluid_model_part.Name)
            variable_utils.SetFlag(flag, value,
                                   self.fluid_model_part.GetSubModelPart(condition_name).Nodes)

    def __CalculateWallDistances(self):
        if (self.distance_calculation_process is None):
            import python_linear_solver_factory as linear_solver_factory
            self.distance_calculation_linear_solver = linear_solver_factory.ConstructSolver(
                self.settings["distance_calculation"]
                ["linear_solver_settings"])
            max_iterations = self.settings["distance_calculation"][
                "max_iterations"].GetInt()
            if (self.domain_size == 2):
                self.distance_calculation_process = Kratos.VariationalDistanceCalculationProcess2D(
                    self.fluid_model_part,
                    self.distance_calculation_linear_solver, max_iterations)
            elif (self.domain_size == 3):
                self.distance_calculation_process = Kratos.VariationalDistanceCalculationProcess3D(
                    self.fluid_model_part,
                    self.distance_calculation_linear_solver, max_iterations)
            else:
                raise Exception("Unsupported domain size")

            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Variational distance calculation process initialized.")

        self.distance_calculation_process.Execute()
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Wall distances calculated.")
