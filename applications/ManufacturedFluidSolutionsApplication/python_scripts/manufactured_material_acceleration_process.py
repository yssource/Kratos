# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS
import KratosMultiphysics.PFEM2Application as PFEM2
import KratosMultiphysics.ConvectionDiffusionApplication as CDA

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
                "use_bfecc"        : false
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.framework = settings["framework"].GetString()
        self.time_scheme = settings["time_scheme"].GetString()
        self.use_bfecc = settings["use_bfecc"].GetBool()

        if (self.time_scheme) and (self.time_scheme != "bdf1"):
            msg = "Requested time scheme: \"" + self.time_scheme + "\""
            msg += "\nAvailable options are:\n"
            msg += "\tNone\n"
            msg += "\t\"bdf1\"\n"
            raise Exception(msg)

        if (self.framework != "eulerian") and (self.framework != "lagrangian"):
            msg = "Requested framework type: \"" + self.framework + "\""
            msg += "\nAvailable options are:\n"
            msg += "\t\"eulerian\"\n"
            msg += "\t\"lagrangian\"\n"
            raise Exception(msg)

        if self.use_bfecc:
            self.locator = KM.BinBasedFastPointLocator2D(self.model_part)
            self.locator.UpdateSearchDatabase()
            self.bfecc = CDA.BFECCConvection2D(self.locator)

    def ExecuteFinalizeSolutionStep(self):
        if self.use_bfecc:
            KM.VariableUtils().CopyVariable(KM.VELOCITY, PFEM2.PROJECTED_VELOCITY, self.model_part.Nodes)
            bfecc_substeps = 10
            self.bfecc.BFECCconvect(self.model_part, PFEM2.PROJECTED_VELOCITY, KM.VELOCITY, bfecc_substeps)
            self.bfecc.CopyVariableToPreviousTimeStep(self.model_part. PFEM2.PROJECTED_VELOCITY)
            KM.VariableUtils().CopyVariable(KM.VELOCITY, PFEM2.PROJECTED_VELOCITY, self.model_part.Nodes)

        if self.time_scheme:
            if self.time_scheme == "bdf1":
                self.manufactured_process.BDF1(KM.VELOCITY, KM.ACCELERATION)

        if self.framework == "eulerian":
            self.manufactured_process.RecoverMaterialAcceleration()            
        if self.framework == "lagrangian":
            KM.VariableUtils().CopyVectorVar(KM.ACCELERATION, KM.MATERIAL_ACCELERATION, self.model_part.Nodes)

        self.manufactured_process.ComputeExactMaterialAcceleration()
        self.manufactured_process.ComputeMaterialAccelerationError()
