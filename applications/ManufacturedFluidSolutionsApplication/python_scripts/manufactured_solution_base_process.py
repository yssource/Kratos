# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS

def Factory(settings, Model):
    raise Exception("Trying the build the ManufacturedSolutionBaseProcess. Please, call the derived classes")

## All the processes python should be derived from "Process"
class ManufacturedSolutionBaseProcess(KM.Process):
    """ Base class for all the sub-processes for the manufactured solutions
    It provides an mandatory method to set the manufactured solution
    and sets the manufactured_process attribute, which is in charge to apply
    the custom velocity, pressure and body force.
    """
    def SetManufacturedSolution(self, manufactured_solution):
        self.manufactured_solution = manufactured_solution

        if hasattr(self, 'model_part'):
            self.manufactured_process = MS.ManufacturedSolutionUtility(self.model_part, self.manufactured_solution)
        else:
            KM.Logger.PrintWarning("ManufacturedSolutionBaseProcess", "No model_part attribute found. The manufactured process is not constructed.")

    def Check(self):
        pass

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def Clear(self):
        pass
