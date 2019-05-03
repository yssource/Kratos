
from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication
try:
    import KratosMultiphysics.MeshMovingApplication
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "not imported")

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_forces_on_nodes_process import ComputeForcesOnNodesProcess
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D

def CreateSolver(cosim_solver_settings, level):
    return KratosPotentialFlowSolver(cosim_solver_settings, level)


class KratosPotentialFlowSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return PotentialFlowAnalysis(self.model, self.project_parameters)

    def InitializeSolutionStep(self):
        super(KratosPotentialFlowSolver, self).InitializeSolutionStep()

        wake_proc_settings = KratosMultiphysics.Parameters(('{'
                                        '"model_part_name": "'
                                        + self.project_parameters["solver_settings"]["ale_boundary_parts"][0].GetString()
                                        + '"}'))
        self.wake_process = DefineWakeProcess2D(self.model, wake_proc_settings)

        conversion_proc_settings = KratosMultiphysics.Parameters(('{'
                                        '"model_part_name": "'
                                        + self.project_parameters["solver_settings"]["ale_boundary_parts"][0].GetString()
                                        + '"}'))
        self.conversion_process = ComputeForcesOnNodesProcess(self.model, conversion_proc_settings)


    def SolveSolutionStep(self):
        self.wake_process.FindWakeElements()
        super(KratosPotentialFlowSolver, self).SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetAnalysisStage().FinalizeSolutionStep()
        self.conversion_process.Execute()
        self.wake_process.CleanMarking()

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__

