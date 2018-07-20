from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication", "FluidDynamicsApplication","ExternalSolversApplication")

# Import applications
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication

# Other imports
import navier_stokes_solver_vmsmonolithic
#import mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return ALENavierStokesSolverVMSMonolithic(model_part, custom_settings)


class ALENavierStokesSolverVMSMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("Monolithic")
        super(ALENavierStokesSolverVMSMonolithic, self).__init__(model_part, navier_stokes_settings)
        # create mesh motion solver
        aux_mesh_solver_settings = KratosMultiphysics.Parameters("""{
            "parallel_type" : "OpenMP"
        }""")
        custom_settings.AddValue("problem_data", aux_mesh_solver_settings)
        print(custom_settings.PrettyPrintJsonString())
        custom_settings.AddValue("solver_settings", custom_settings["ale_settings"])
        custom_settings.RemoveValue("ale_settings")
        # err

        ## TODO add warning if "compute_reactions" is set to false!

        import python_solvers_wrapper_mesh_motion
        self.ale_solver = python_solvers_wrapper_mesh_motion.CreateSolver(self.main_model_part, custom_settings)

        print("::[ALENavierStokesSolverVMSMonolithic]:: Construction finished")

    def AddVariables(self):
        super(ALENavierStokesSolverVMSMonolithic, self).AddVariables()
        self.ale_solver.AddVariables()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Variables ADDED.")

    def AddDofs(self):
        super(ALENavierStokesSolverVMSMonolithic, self).AddDofs()
        self.ale_solver.AddDofs()
        print("::[ALENavierStokesSolverVMSMonolithic]:: DOFs ADDED.")

    def Initialize(self):
        super(ALENavierStokesSolverVMSMonolithic, self).Initialize()
        self.ale_solver.Initialize()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALENavierStokesSolverVMSMonolithic, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def SolveSolutionStep(self):
        self.GetMeshMotionSolver().SolveSolutionStep()
        # Copy the MESH_VELOCITY to the VELOCITY (ALE) on the interface
        for i in range(self.ale_solver.settings["ale_interface_parts"].size()):
            part_name = self.ale_solver.settings["ale_interface_parts"][i].GetString()
            part_nodes = self.main_model_part.GetSubModelPart(part_name).Nodes
            KratosMultiphysics.VariableUtils().CopyVectorVar(KratosMultiphysics.MESH_VELOCITY, KratosMultiphysics.VELOCITY, part_nodes)

        self.GetFluidSolver().SolveSolutionStep()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
