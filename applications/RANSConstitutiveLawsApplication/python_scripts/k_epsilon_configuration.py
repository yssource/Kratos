from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.RANSConstitutiveLawsApplication as RANS

from KratosMultiphysics.kratos_utilities import IsApplicationAvailable
if IsApplicationAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication
else:
    raise Exception("k_epsilon_configuration.py depends on the FluidDynamicsApplication, which is not available.")

def Factory(model, parameters):


class KEpsilonConfiguration(FluidDynamicsApplication.TurbulenceModelConfiguration):

    def __init__(self, model, parameters):
        super(KEpsilonConfiguration,self).__init__(self,model,parameters)
        self.default_settings = KratosMultiphysics.Parameters(r'''{
            "model_type"            : "k_epsilon",
            "inlet_conditions"      : [],
            "outlet_conditions"     : [],
            "wall_conditions"       : [],
            "max_distance_calculation_iterations" : 10,
            "mesh_moving"       : false,
            "echo_level"        : 0,
            "model_properties"  : {
                "time_scheme"    : "transient",
                "scheme_settings": {
                    "scheme_type": "bossak",
                    "alpha_bossak": -0.3
                },
                "echo_level":0,
                "show_warnings"     : false,
                "convergence_tolerances":
                {
                    "k_relative_tolerance": 1e-3,
                    "k_absolute_tolerance": 1e-5,
                    "epsilon_relative_tolerance": 1e-3,
                    "epsilon_absolute_tolerance": 1e-5,
                    "turbulent_viscosity_relative_tolerance": 1e-3,
                    "turbulent_viscosity_absolute_tolerance": 1e-5,
                    "k_max_iterations": 200,
                    "epsilon_max_iterations": 200,
                    "maximum_coupling_iterations": 200,
		            "maximum_stabilization_multiplier": 1e+15,
                    "echo_level": 0
                },
                "flow_parameters":
                {
                    "ramp_up_time"                : 0.5,
                    "turbulent_mixing_length"     : 0.07,
                    "turbulence_intensity"        : 0.055175811859971446,
                    "turbulent_viscosity_min"     : 1e-12,
                    "turbulent_viscosity_max"     : 1e+2
                }
            }
        }''')

        self.settings = self.ValidateInputSettings(parameters, default_settings)

        self.rans_elements = ["RANSEVMK","RANSEVMEPSILON"]
        self.rans_conditions = ["Condition","Condition"]


    def AddVariables(self):
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2)
        self.fluid_model_part.AddNodalSolutionStepVariable(RANS_Y_PLUS)
        self.fluid_model_part.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1)
        self.fluid_model_part.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2)

    def AddDofs(self):
        self.fluid_model_part.AddDofs(TURBULENT_KINETIC_ENERGY)
        self.fluid_model_part.AddDofs(TURBULENT_ENERGY_DISSIPATION_RATE)

    def GetTurbulenceModelProcess(self):
        k_solver = None
        epsilon_solver = None
        distance_solver = None
        if self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            return RANS.TurbulenceEvmKEpsilon2DProcess(self.settings, k_solver, epsilon_solver, distance_solver)
        else:
            return RANS.TurbulenceEvmKEpsilon3DProcess(self.settings, k_solver, epsilon_solver, distance_solver)
