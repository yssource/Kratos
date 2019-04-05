# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS

from KratosMultiphysics.ManufacturedFluidSolutionsApplication.manufactured_solution_base_process import ManufacturedSolutionBaseProcess as ManufacturedProcess

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomVelocityProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomVelocityProcess(ManufacturedProcess):
    '''
    This process apply the custom velocity to the conditions.
    The velocity field is defined by the manufactured solution.
    TODO: use the parameters previously stored in the ProceesInfo
    '''
    def __init__(self, model, settings ):

        default_settings = KM.Parameters("""
            {
                "model_part_name"          : "model_part_name",
                "constrained"              : [true, true, true]
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # The model part manufactured solution applies to
        self.model_part = model[settings["model_part_name"].GetString()]

        # Fixity process
        self.is_fixed = [False]*3
        self.is_fixed[0] = settings["constrained"][0].GetBool()
        self.is_fixed[1] = settings["constrained"][1].GetBool()
        self.is_fixed[2] = settings["constrained"][2].GetBool()
        self.fix_variable = [None]*3
        self.fix_variable[0] = KM.VELOCITY_X
        self.fix_variable[1] = KM.VELOCITY_Y
        self.fix_variable[2] = KM.VELOCITY_Z
        self.variable_utils = KM.VariableUtils()


    def ExecuteBeforeSolutionLoop(self):
        self.manufactured_process.SetVelocity()
        for i in range(3):
            self.variable_utils.ApplyFixity(self.fix_variable[i], self.is_fixed[i], self.model_part.Nodes)


    def ExecuteInitializeSolutionStep(self):
        self.manufactured_process.SetVelocity()
        for i in range(3):
            self.variable_utils.ApplyFixity(self.fix_variable[i], self.is_fixed[i], self.model_part.Nodes)
