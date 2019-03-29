# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionApplication as MS


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomBodyForceProcess(KM.Process):
    def __init__(self, model, settings ):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"          : "model_part_name",
                "manufactured_name"        : "manufactured_solution_name",
                "manufactured_parameters"  : {},
                "set_initial_values"       : True,
                "compute_relative_error"   : True
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # We construct the manufactured solution
        fluid_property = self.model_part.ElementsArray(0)[0].GetProperties()
        manufactured_class = getattr(MS, settings["manufactured_name"].GetString())
        self.manufactured = manufactured_class(fluid_properties, settings["manufactured_parameters"])

        # We construct the process to apply the manufactured solution
        self.model_part = model[settings["model_part_name"].GetString()]
        self.manufactured_process = MS.ManufacturedSolutionUtility(self.model_part, self.manufactured)

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
