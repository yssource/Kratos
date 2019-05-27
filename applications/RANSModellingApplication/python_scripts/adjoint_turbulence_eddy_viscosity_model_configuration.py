import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSModellingApplication as KratosRANS
import math

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
    from adjoint_turbulence_model_configuration import AdjointTurbulenceModelConfiguration
else:
    msg = "RANSModellingApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)


class AdjointTurbulenceEddyViscosityModelConfiguration(
        AdjointTurbulenceModelConfiguration):
    def __init__(self, model, settings):
        self._validate_settings_in_baseclass = True  # To be removed eventually

        super(AdjointTurbulenceEddyViscosityModelConfiguration, self).__init__(
            model, settings)

        # self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.adjoint_y_plus_model_process = None
        self.element_name = None

    def GetDefaultSettings(self):
        return Kratos.Parameters(r'''{
            "model_type"            : "",
            "model_settings"        : {},
            "adjoint_y_plus_model"  : {
                "model_type"     : "",
                "model_settings" : {}
            },
            "echo_level"              : 0
        }''')

    def Initialize(self):
        super(AdjointTurbulenceEddyViscosityModelConfiguration,
              self).Initialize()

        self.PrepareSolvingStrategy()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Initialization successfull.")

    def AddVariables(self):
        # adding variables required by rans eddy viscosity models
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.KINEMATIC_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.TURBULENT_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.RANS_Y_PLUS)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.RELAXED_ACCELERATION)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Successfully added solution step variables.")

    def GetAdjointYPlusModel(self):
        if self.adjoint_y_plus_model_process is None:
            y_plus_model_settings = self.settings["adjoint_y_plus_model"]
            import adjoint_rans_y_plus_model_factory
            adjoint_rans_y_plus_model_factory.InitializeModelPartName(
                y_plus_model_settings, self.model, self.fluid_model_part)
            self.adjoint_y_plus_model_process = adjoint_rans_y_plus_model_factory.Factory(
                y_plus_model_settings, self.model)
            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Initialized " + self.adjoint_y_plus_model_process.__str__())

        return self.adjoint_y_plus_model_process

    def GetAdjointResponseFunction(self, settings):
        domain_size = self.fluid_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

        if (settings["response_type"].GetString() == "drag"):
            if (domain_size == 2):
                return KratosRANS.RansDragResponseFunction2D(settings["custom_settings"], self.fluid_model_part)
            elif (domain_size == 3):
                return KratosRANS.RansDragResponseFunction3D(settings["custom_settings"], self.fluid_model_part)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("Unknown RANS response_type \"" + settings["response_type"].GetString() + "\"")
