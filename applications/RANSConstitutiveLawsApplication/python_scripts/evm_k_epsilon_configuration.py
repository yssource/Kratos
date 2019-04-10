from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.RANSConstitutiveLawsApplication as RANS
from turbulence_eddy_viscosity_model_configuration import TurbulenceEddyViscosityModelConfiguration

from KratosMultiphysics.kratos_utilities import IsApplicationAvailable
if IsApplicationAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication
else:
    raise Exception("k_epsilon_configuration.py depends on the FluidDynamicsApplication, which is not available.")



class TurbulenceKEpsilonConfiguration(TurbulenceEddyViscosityModelConfiguration):

    def __init__(self, model, parameters):
        super(TurbulenceKEpsilonConfiguration,self).__init__(self,model,parameters)

        default_settings = KM.Parameters(r'''{
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
                "coupling_settings":{
                    "relative_tolerance"    : 1e-3,
                    "absolute_tolerance"    : 1e-5,
                    "max_iterations"        : 200,
                    "echo_level"            : 0                    
                }          
                "flow_parameters":
                {
                    "ramp_up_time"                        : 0.5,
                    "turbulent_mixing_length"             : 0.07,
                    "turbulence_intensity"                : 0.055175811859971446,
                    "turbulent_viscosity_min_max_factor"  : 1e-6
                }
        }''')

        parameters["model_settings"].ValidateAndAssignDefaults(default_settings)  

        self.settings = self.ValidateInputSettings(parameters, default_settings)

        self.rans_elements = ["RANSEVMK","RANSEVMEPSILON"]
        self.rans_conditions = ["Condition","Condition"]
        self.rans_solvers = []
        self.linear_solvers = []
        self.conv_criterians = []

    def PrepareSolvers(self):




    def AddVariables(self):
        super(TurbulenceKEpsilonConfiguration, self).AddVariables()

        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE)
        self.fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2)
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
