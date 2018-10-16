from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")
import KratosMultiphysics.MeshMovingApplication

# Other imports
from ale_fluid_solver import ALEFluidSolver
import python_solvers_wrapper_mesh_motion


def CreateSolver(model, solver_settings, parallelism):
    return ALEChimeraSolver(model, solver_settings, parallelism)

class ALEChimeraSolver(ALEFluidSolver):
    def __init__(self,model, solver_settings, parallelism):
        super(ALEChimeraSolver,self).__init__(model, solver_settings, parallelism)

    def _CreateFluidSolver(self, solver_settings, parallelism):
        KratosMultiphysics.CheckRegisteredApplications("ChimeraApplication")
        import python_solvers_wrapper_fluid_chimera
        return python_solvers_wrapper_fluid_chimera.CreateSolverByParameters(
            self.model, solver_settings, parallelism)
        