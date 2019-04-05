# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS

from manufactured_solution_base_process import ManufacturedSolutionBaseProcess as ManufacturedProcess

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomPressureProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomPressureProcess(ManufacturedProcess):
    '''
    This process apply the custom velocity to the conditions.
    The velocity field is defined by the manufactured solution.
    TODO: use the parameters previously stored in the ProceesInfo
    '''
    def __init__(self, model, settings ):

        default_settings = KM.Parameters("""
            {
                "model_part_name"          : "model_part_name",
                "constrained"              : true,
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # The model part manufactured solution applies to
        self.model_part = model[settings["model_part_name"].GetString()]

        # Fixity process
        self.is_fixed = settings["constrained"].GetBool()
        self.variable_utils = KM.VariableUtils()


    def ExecuteBeforeSolutionLoop(self):
        self.manufactured_process.SetPressure()


    def ExecuteInitializeSolutionStep(self):
        self.manufactured_process.SetPressure()
        self.variable_utils.ApplyFixity(KM.PRESSURE, self.is_fixed, self.model_part.Nodes)
