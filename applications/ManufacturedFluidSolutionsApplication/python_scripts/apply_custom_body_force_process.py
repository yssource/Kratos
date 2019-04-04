# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionApplication as MS

from manufactured_solution_base_process import ManufacturedSolutionBaseProcess as ManufacturedProcess

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomBodyForceProcess(Model, settings["Parameters"])

## All the python manufactured processes should be derived from a base class, which is derived from "Process"
class ApplyCustomBodyForceProcess(ManufacturedProcess):
    '''
    This process apply the custom body force to all the fluid domain
    The body force is defined by the manufactured solution
    TODO: use the parameters previously stored in the ProceesInfo
    '''
    def __init__(self, model, settings ):

        default_settings = KM.Parameters("""
            {
                "model_part_name"          : "model_part_name",
                "set_initial_values"       : True,
                "compute_relative_error"   : True
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # The model part manufactured solution applies to
        self.model_part = model[settings["model_part_name"].GetString()]

        # Auxiliary variables
        self.set_initial_values = settings["set_initial_values"].GetBool()
        self.compute_error = settings["copmute_relative_error"].GetBool()


    def ExecuteBeforeSolutionLoop(self):
        if self.set_initial_values:
            self.manufactured_process.SetVelocity()
            self.manufactured_process.SetPressure()
        if self.compute_error:
            self.manufactured_process.ComputeExactVelocity()
            self.manufactured_process.ComputeExactPressure()
            self.manufactured_procees.ComputeVelocityRelativeError()
            self.manufactured_process.ComputePressureRelativeError()


    def ExecuteInitializeSolutionStep(self):
        self.manufactured_process.SetBodyForce()


    def ExecuteBeforeOutputStep(self):
        if self.compute_error:
            self.manufactured_process.ComputeExactVelocity()
            self.manufactured_process.ComputeExactPressure()
            self.manufactured_procees.ComputeVelocityRelativeError()
            self.manufactured_process.ComputePressureRelativeError()
