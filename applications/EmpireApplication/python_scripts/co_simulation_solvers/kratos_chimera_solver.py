from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ChimeraApplication
try:
    import KratosMultiphysics.MeshMovingApplication
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "not imported")

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from fluid_dynamics_analysis import FluidDynamicsAnalysis
from fluid_chimera_analysis import FluidChimeraAnalysis

def CreateSolver(cosim_solver_settings, level):
    return KratosFluidSolver(cosim_solver_settings, level)

class KratosFluidSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return FluidChimeraAnalysis(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()
    
    def SolveSolutionStep(self):
        self._GetAnalysisStage()._GetSolver().SolveSolutionStep()

    def _Name(self):
        return self.__class__.__name__
