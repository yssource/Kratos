# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS
# import KratosMultiphysics.DEMApplication as DEM
# import KratosMultiphysics.SwimmingDEMApplication as SD

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
                "framework"        : "eulerian",
                "time_scheme"      : ""
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.framework = settings["framework"].GetString()
        self.time_scheme = settings["time_scheme"].GetString()

        if self.time_scheme and self.time_scheme is not "bdf1":
            msg = "Requested time scheme: " + self.time_scheme
            msg += "\nAvailable options are:\n"
            msg += "\tNone\n"
            msg += "\t\"bdf1\"\n"

        if self.framework == "eulerian":
            # recovery_parameters = KM.Parameters()
            # sub_param = recovery_parameters.AddEmptyValue("variables_for_recovery")
            # sub_param = sub_param.AddEmptyValue("material_derivative")
            # sub_param = sub_param.AddEmptyValue("VELOCITY").SetString("MATERIAL_ACCELERATION")
            # self.derivative_recoverer = SD.StandardRecoveryUtility(self.model_part, recovery_parameters)
            # KM.CalculateNodalAreaProcess(self.model_part)
            pass
        elif self.framework == "lagrangian":
            pass
        else:
            msg = "Requested framework type: " + self.framework
            msg += "\nAvailable options are:\n"
            msg += "\t\"eulerian\"\n"
            msg += "\t\"lagrangian\"\n"
            raise Exception(msg)

    def ExecuteFinalizeSolutionStep(self):
        if self.time_scheme:
            if self.time_scheme == "bdf1":
                self.manufactured_process.BDF1(KM.VELOCITY, KM.ACCELERATION)
        if self.framework == "eulerian":
            # self.derivative_recoverer.Recover()
            self.manufactured_process.RecoverMaterialAcceleration()            
        if self.framework == "lagrangian":
            KM.VariableUtils().CopyVectorVar(KM.ACCELERATION, KM.MATERIAL_ACCELERATION, self.model_part.Nodes)
        self.manufactured_process.ComputeExactMaterialAcceleration()
        self.manufactured_process.ComputeMaterialAccelerationError()
