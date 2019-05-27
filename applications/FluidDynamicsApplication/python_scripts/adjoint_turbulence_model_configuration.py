from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
from python_solver import PythonSolver
from kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("RANSModellingApplication"):
    import KratosMultiphysics.RANSModellingApplication as KratosRANS


def CreateTurbulenceModel(model, settings):
    if not CheckIfApplicationsAvailable("RANSModellingApplication"):
        msg = "Using an adjoint turbulence model requires the RANSModellingApplication. "
        msg += "Please re-install/re-compile with RANSModellingApplication."
        raise Exception(msg)

    from adjoint_turbulence_model_factory import Factory
    return Factory(settings, model)

class AdjointTurbulenceModelConfiguration(PythonSolver):
    def AddVariables(self):
        msg = "Calling the base AdjointTurbulenceModelConfiguration class AddVariables method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def AddDofs(self):
        msg = "Calling the base AdjointTurbulenceModelConfiguration class AddDofs method."
        msg += " Please override it in the derrived class."
        raise Exception(msg)

    def GetAdjointElementName(self):
        msg = "Calling the base AdjointTurbulenceModelConfiguration class GetAdjointElementName method."
        msg += " Please override it in the derrived class to return an element name."
        raise Exception(msg)

    def GetAdjointResponseFunction(self, settings):
        msg = "Calling the base AdjointTurbulenceModelConfiguration class GetAdjointElementName method."
        msg += " Please override it in the derrived class to return a response function."
        raise Exception(msg)

    def GetAdjointYPlusModel(self):
        msg = "Calling the base AdjointTurbulenceModelConfiguration class GetAdjointYPlusModel method."
        msg += " Please override it in the derrived class to return a process."
        raise Exception(msg)

    def Initialize(self):
        self.GetAdjointYPlusModel().ExecuteInitialize()

    def Check(self):
        self.GetAdjointYPlusModel().Check()

    def InitializeSolutionStep(self):
        self.GetAdjointYPlusModel().ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.GetAdjointYPlusModel().ExecuteFinalizeSolutionStep()