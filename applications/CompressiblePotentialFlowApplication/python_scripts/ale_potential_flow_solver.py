from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.MeshMovingApplication.ale_fluid_solver import AleFluidSolver
import KratosMultiphysics.MeshMovingApplication.python_solvers_wrapper_mesh_motion as mesh_mothion_solvers_wrapper
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver as potential_flow_solver

def CreateSolver(model, solver_settings, parallelism):
    return AlePotentialFlowSolver(model, solver_settings, parallelism)


class AlePotentialFlowSolver(AleFluidSolver):
    def __init__(self, model, solver_settings, parallelism):
        default_settings = KratosMultiphysics.Parameters("""{
            "solver_type"                 : "ale_fluid",
            "echo_level"                  : 0,
            "ale_boundary_parts"          : [ ],
            "mesh_motion_parts"           : [ ],
            "fluid_solver_settings"       : { },
            "mesh_motion_solver_settings" : { }
        }""")
        solver_settings.ValidateAndAssignDefaults(default_settings)
        self.parallelism = parallelism

        super(AleFluidSolver, self).__init__(model, solver_settings) # TODO: Isn't there a better way to call grandparent init ?

        ## Creating the fluid part
        fluid_solver_settings = self.settings["fluid_solver_settings"]
        fluid_model_part_name = fluid_solver_settings["model_part_name"].GetString()
        if not self.model.HasModelPart(fluid_model_part_name):
            model.CreateModelPart(fluid_model_part_name)
        self.fluid_solver = self._CreateFluidSolver(fluid_solver_settings, parallelism)

        [fluid_solver_settings, mesh_motion_solver_settings] = self._CheckSettingsConsistency(fluid_solver_settings, self.settings["mesh_motion_solver_settings"])

        ## Creating the mesh-motion solver
        if not mesh_motion_solver_settings.Has("echo_level"):
            mesh_motion_solver_settings.AddValue("echo_level", self.settings["echo_level"])
        self.mesh_motion_solver_full_mesh = mesh_mothion_solvers_wrapper.CreateSolverByParameters(
            model, mesh_motion_solver_settings, parallelism)

        self.fluid_solver.min_buffer_size = 2
        KratosMultiphysics.Logger.PrintInfo("::[AlePotentialFlowSolver]::", "Construction finished")

    def AddVariables(self):
        self.mesh_motion_solver_full_mesh.AddVariables()
        self.fluid_solver.AddVariables()
        KratosMultiphysics.Logger.PrintInfo("::[AlePotentialFlowSolver]::", "Variables Added")

    def _CreateFluidSolver(self, solver_settings, parallelism):
        return potential_flow_solver.CreateSolver(self.model, solver_settings)

    def AdvanceInTime(self, current_time):
        return 0.0 # potential flow is steady state

    def SolveSolutionStep(self):
        for mesh_solver in self.mesh_motion_solvers:
            mesh_solver.SolveSolutionStep()
        self.fluid_solver.SolveSolutionStep()