# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionApplication as MS


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomPressureProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomPressureProcess(KM.Process):
    '''
    This process apply the custom velocity to the conditions.
    The velocity field is defined by the manufactured solution.
    TODO: use the parameters previously stored in the ProceesInfo
    '''
    def __init__(self, model, settings ):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"          : "model_part_name",
                "manufactured_name"        : "manufactured_solution_name",
                "manufactured_parameters"  : {},
                "constrained"              : true,
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # We construct the manufactured solution
        self.model_part = model[settings["model_part_name"].GetString()]
        fluid_property = self.model_part.ElementsArray(0)[0].GetProperties()
        manufactured_class = getattr(MS, settings["manufactured_name"].GetString())
        self.manufactured = manufactured_class(fluid_properties, settings["manufactured_parameters"])

        # We construct the process to apply the manufactured solution
        self.manufactured_process = MS.ManufacturedSolutionUtility(self.model_part, self.manufactured)

        # Fixity process
        self.is_fixed = settings["constrained"].GetBool()
        self.variable_utils = KM.VariableUtils()


    def ExecuteBeforeSolutionLoop(self):
        self.manufactured_process.SetPressure()


    def ExecuteInitializeSolutionStep(self):
        self.manufactured_process.SetPressure()
        self.variable_utils.ApplyFixity(KM.PRESSURE, self.is_fixed, self.model_part.Nodes)
