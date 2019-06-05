# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SD

from KratosMultiphysics.ManufacturedFluidSolutionsApplication.manufactured_solution_base_process import ManufacturedSolutionBaseProcess as ManufacturedBaseProcess

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ManufacturedMaterialAccelerationProcess(Model, settings["Parameters"])

## All the python manufactured processes should be derived from a base class, which is derived from "Process"
class ManufacturedMaterialAccelerationProcess(ManufacturedBaseProcess):
    def __init__(self, model, settings ):

        default_settings = KM.Parameters("""
            {
                "model_part_name"  : "model_part_name",
                "framework"        : "eulerian"
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # The model part manufactured solution applies to
        self.model_part = model[settings["model_part_name"].GetString()]

        self.framework = settings["framework"].GetString()
    
        if self.framework == "eulerian":
            recovery_parameters = KM.Parameters()
            sub_param = recovery_parameters.AddEmptyValue("variables_for_recovery")
            sub_param = sub_param.AddEmptyValue("material_derivative")
            sub_param = sub_param.AddEmptyValue("VELOCITY").SetString("MATERIAL_ACCELERATION")
            self.derivative_recoverer = SD.StandardRecoveryUtility(self.model_part, recovery_parameters)
        if self.framework == "lagrangian":
            self.velocity_variable = KM.KratosGlobals.GetVariable("VELOCITY")
            self.acceleration_variable = KM.KratosGlobals.GetVariable("MATERIAL_ACCELERATION")

    def ExecuteFinalizeSolutionStep(self):
        if self.framework == "eulerian":
            self.derivative_recoverer.Recover()
        if self.framework == "lagrangian":
            KM.VariableUtils().CopyVectorVar(self.velocity_variable, self.acceleration_variable, self.model_part.Nodes)
        self.manufactured_process.CopmuteMaterialAccelerationError
