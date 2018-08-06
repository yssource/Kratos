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


def CreateSolver(model, custom_settings):
    return ALENavierStokesSolverVMSMonolithic(model, custom_settings)


class ALENavierStokesSolverVMSMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic):
    def __init__(self, model, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        fluid_settings = custom_settings.Clone()
        fluid_solver_settings = fluid_settings["solver_settings"]
        fluid_solver_settings.RemoveValue("ale_settings")

        mesh_motion_settings = custom_settings.Clone()
        mesh_motion_settings.RemoveValue("solver_settings")
        mesh_motion_settings.AddValue("solver_settings", custom_settings["solver_settings"]["ale_settings"])
        mesh_motion_solver_settings = mesh_motion_settings["solver_settings"]

        fluid_model_part_name = fluid_solver_settings["model_part_name"].GetString()
        if not model.HasModelPart(fluid_model_part_name):
            fluid_mesh_model_part = KratosMultiphysics.ModelPart(fluid_model_part_name)
            model.AddModelPart(fluid_mesh_model_part)

        super(ALENavierStokesSolverVMSMonolithic, self).__init__(model, fluid_solver_settings)

        # create mesh motion solver
        if not mesh_motion_solver_settings.Has("echo_level"):
            mesh_motion_solver_settings.AddEmptyValue("echo_level")
            fluid_echo_lvl = fluid_solver_settings["echo_level"].GetInt()
            mesh_motion_solver_settings["echo_level"].SetInt(fluid_echo_lvl)

        if mesh_motion_solver_settings.Has("model_part_name"):
            if not fluid_model_part_name == mesh_motion_solver_settings["model_part_name"].GetString():
                raise Exception('Fluid- and Mesh-Solver have to use the same "model_part_name"!')
        else:
            mesh_motion_solver_settings.AddEmptyValue("model_part_name")
            fluid_model_part_name = fluid_solver_settings["model_part_name"].GetString()
            mesh_motion_solver_settings["model_part_name"].SetString(fluid_model_part_name)

        fluid_domain_size = fluid_solver_settings["domain_size"].GetInt()
        if mesh_motion_solver_settings.Has("domain_size"):
            mesh_motion_domain_size = mesh_motion_solver_settings["domain_size"].GetInt()
            if not fluid_domain_size == mesh_motion_domain_size:
                raise Exception('Fluid- and Mesh-Solver have to use the same "domain_size"!')
        else:
            mesh_motion_solver_settings.AddEmptyValue("domain_size")
            mesh_motion_solver_settings["domain_size"].SetInt(fluid_domain_size)

        self.ale_interface_part_names = []
        if mesh_motion_solver_settings.Has("ale_interface_parts"):
            for i in range(mesh_motion_solver_settings["ale_interface_parts"].size()):
                self.ale_interface_part_names.append(
                    mesh_motion_solver_settings["ale_interface_parts"][i].GetString())
            mesh_motion_solver_settings.RemoveValue("ale_interface_parts")

        import python_solvers_wrapper_mesh_motion
        self.mesh_motion_solver = python_solvers_wrapper_mesh_motion.CreateSolver(model, mesh_motion_settings)

        # Getting the min_buffer_size from both solvers
        # and assigning it to the fluid_solver, bcs this one handles the model_part
        # self.min_buffer_size = max(super(ALENavierStokesSolverVMSMonolithic, self).GetMinimumBufferSize(),
        #                                         self.mesh_motion_solver.GetMinimumBufferSize())

        self.is_printing_rank = super(ALENavierStokesSolverVMSMonolithic, self)._IsPrintingRank()

        # TODO move to "Check"?
        if (self.mesh_motion_solver.settings["calculate_mesh_velocities"].GetBool() == False
            and self.is_printing_rank):
            info_msg = "Mesh velocities are not being computed in the Mesh solver!"
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        if self.settings["compute_reactions"].GetBool() == False and self.is_printing_rank:
            info_msg = "Reactions are not being computed in the Fluid solver!"
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", info_msg)

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Construction finished")

    def AddVariables(self):
        self.mesh_motion_solver.AddVariables()
        super(ALENavierStokesSolverVMSMonolithic, self).AddVariables()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Variables ADDED.")

    def AddDofs(self):
        self.mesh_motion_solver.AddDofs()
        super(ALENavierStokesSolverVMSMonolithic, self).AddDofs()
        print("::[ALENavierStokesSolverVMSMonolithic]:: DOFs ADDED.")

    def ImportModelPart(self):
        super(ALENavierStokesSolverVMSMonolithic, self).ImportModelPart() # only ONE solver imports the ModelPart

    def PrepareModelPart(self):
        # Doing it ONLY for the fluid solver (since this contains filling the buffer)
        super(ALENavierStokesSolverVMSMonolithic, self).PrepareModelPart()

    def AdvanceInTime(self, current_time):
        # Doing it ONLY for the fluid solver
        return super(ALENavierStokesSolverVMSMonolithic, self).AdvanceInTime(current_time) # returning new time

    def Initialize(self):
        self.mesh_motion_solver.Initialize()
        super(ALENavierStokesSolverVMSMonolithic, self).Initialize()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Finished initialization.")

    def Finalize(self):
        self.mesh_motion_solver.Finalize()
        super(ALENavierStokesSolverVMSMonolithic, self).Finalize()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Finished initialization.")

    def InitializeSolutionStep(self):
        self.mesh_motion_solver.InitializeSolutionStep()
        super(ALENavierStokesSolverVMSMonolithic, self).InitializeSolutionStep()

    def Predict(self):
        self.mesh_motion_solver.Predict()
        super(ALENavierStokesSolverVMSMonolithic, self).Predict()

    def FinalizeSolutionStep(self):
        self.mesh_motion_solver.FinalizeSolutionStep()
        super(ALENavierStokesSolverVMSMonolithic, self).FinalizeSolutionStep()

    def SolveSolutionStep(self):
        self.GetMeshMotionSolver().SolveSolutionStep()
        # Copy the MESH_VELOCITY to the VELOCITY (ALE) on the interface
        for part_name in self.ale_interface_part_names:
            part_nodes = self.GetComputingModelPart().GetSubModelPart(part_name).Nodes
            KratosMultiphysics.VariableUtils().CopyVectorVar(KratosMultiphysics.MESH_VELOCITY,
                                                             KratosMultiphysics.VELOCITY,
                                                             part_nodes)

        self.GetFluidSolver().SolveSolutionStep()

    def GetFluidSolver(self):
        return super(ALENavierStokesSolverVMSMonolithic, self)

    def GetMeshMotionSolver(self):
        return self.mesh_motion_solver

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()

    def GetComputingModelPart(self):
        return self.mesh_motion_solver.GetComputingModelPart() # this is the same as the one used in Fluid

    def Check(self):
        self.mesh_motion_solver.Check()
        super(ALENavierStokesSolverVMSMonolithic, self).Check()

    def Clear(self):
        self.mesh_motion_solver.Clear()
        super(ALENavierStokesSolverVMSMonolithic, self).Clear()