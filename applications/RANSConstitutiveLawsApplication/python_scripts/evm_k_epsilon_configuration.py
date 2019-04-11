from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSConstitutiveLawsApplication as KratosRANS
from turbulence_eddy_viscosity_model_configuration import TurbulenceEddyViscosityModelConfiguration

from KratosMultiphysics.kratos_utilities import IsApplicationAvailable
if IsApplicationAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication
else:
    raise Exception(
        "k_epsilon_configuration.py depends on the FluidDynamicsApplication, which is not available."
    )

class TurbulenceKEpsilonConfiguration(
        TurbulenceEddyViscosityModelConfiguration):
    def __init__(self, model, parameters):
        super(TurbulenceKEpsilonConfiguration, self).__init__(model, parameters)

        default_settings = Kratos.Parameters(r'''{
                "time_scheme"    : "transient",
                "scheme_settings": {
                    "scheme_type": "bossak",
                    "alpha_bossak": -0.3
                },
                "echo_level"        :0,
                "show_warnings"     : false,
                "turbulent_kinetic_energy_settings":{
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "echo_level"            : 0,
                    "linear_solver_settings": {}
                },
                "turbulent_energy_dissipation_rate_settings":{
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "echo_level"            : 0,
                    "linear_solver_settings": {}
                },
                "flow_parameters":
                {
                    "ramp_up_time"                        : 0.5,
                    "turbulent_mixing_length"             : 0.07,
                    "turbulence_intensity"                : 0.055175811859971446,
                    "turbulent_viscosity_min_max_factor"  : 1e-6
                }
        }''')

        parameters["model_settings"].ValidateAndAssignDefaults(
            default_settings)

        self.model_elements = ["RANSEVMK", "RANSEVMEPSILON"]
        self.model_conditions = ["Condition", "Condition"]
        self.rans_solver_configurations = []
        self.is_initial_values_assigned = False

        self.ramp_up_time = self.settings["model_settings"]["flow_parameters"]["ramp_up_time"].GetDouble()

    def PrepareSolvingStrategy(self):
        scheme_settings = self.settings["model_settings"]["scheme_settings"]

        # create turbulent kinetic energy strategy
        model_part = self.turbulence_model_parts_list[0]
        solver_settings = self.settings["model_settings"]["turbulent_kinetic_energy_settings"]
        scalar_variable = KratosRANS.TURBULENT_KINETIC_ENERGY
        scalar_variable_rate = KratosRANS.TURBULENT_KINETIC_ENERGY_RATE
        relaxed_scalar_variable_rate = KratosRANS.RANS_AUXILIARY_VARIABLE_1
        self.rans_solver_configurations.append(
            self.CreateStrategy(solver_settings, scheme_settings, model_part,
                                scalar_variable, scalar_variable_rate,
                                relaxed_scalar_variable_rate))
        self.strategies_list.append(self.rans_solver_configurations[-1][0])

        # create turbulent energy dissipation rate strategy
        model_part = self.turbulence_model_parts_list[1]
        solver_settings = self.settings["model_settings"]["turbulent_energy_dissipation_rate_settings"]
        scalar_variable = KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE
        scalar_variable_rate = KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2
        relaxed_scalar_variable_rate = KratosRANS.RANS_AUXILIARY_VARIABLE_2
        self.rans_solver_configurations.append(
            self.CreateStrategy(solver_settings, scheme_settings, model_part,
                                scalar_variable, scalar_variable_rate,
                                relaxed_scalar_variable_rate))
        self.strategies_list.append(self.rans_solver_configurations[-1][0])

    def AddVariables(self, model_part):
        super(TurbulenceKEpsilonConfiguration, self).AddVariables(model_part)

        # adding k-epsilon specific variables
        model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2)
        model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)

    def AddDofs(self, model_part):
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, model_part)
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, model_part)

    def InitializeSolutionStep(self):
        time = self.fluid_model_part.GetProcessInfo()[Kratos.TIME]
        if (time >= self.ramp_up_time):
            self.is_computing_solution = True
            super(TurbulenceKEpsilonConfiguration, self).InitializeSolutionStep()

    def InitializeBoundaryNodes(self):
        rans_variable_utils = KratosRANS.RansVariableUtils()

        rans_variable_utils.FixScalarVariableDofs(Kratos.INLET, KratosRANS.TURBULENT_KINETIC_ENERGY, self.fluid_model_part)
        rans_variable_utils.FixScalarVariableDofs(Kratos.INLET, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.fluid_model_part)

        rans_variable_utils.FixScalarVariableDofs(Kratos.STRUCTURE, KratosRANS.TURBULENT_KINETIC_ENERGY, self.fluid_model_part)
        rans_variable_utils.FixScalarVariableDofs(Kratos.STRUCTURE, KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.fluid_model_part)

    def UpdateTurbulentViscosity(self):
        evm_k_epsilon_utilities = KratosRANS.EvmKepsilonModelUtilities()
        evm_k_epsilon_utilities.CalculateTurbulentViscosity(self.fluid_model_part)

    def UpdateBoundaryConditions(self):
        self.CalculateYPlus()
        evm_k_epsilon_utilities = KratosRANS.EvmKepsilonModelUtilities()
        evm_k_epsilon_utilities.UpdateBoundaryConditions(self.fluid_model_part)

        if not self.is_initial_values_assigned:
            evm_k_epsilon_utilities.AssignInitialValues(self.fluid_model_part)
            self.is_initial_values_assigned = True
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Assigned initial values for turbulence model.")

