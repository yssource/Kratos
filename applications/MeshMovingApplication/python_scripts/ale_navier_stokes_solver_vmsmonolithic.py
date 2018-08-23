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
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")

        fluid_model_part_name = navier_stokes_settings["model_part_name"].GetString()
        if not model.HasModelPart(fluid_model_part_name):
            fluid_mesh_model_part = KratosMultiphysics.ModelPart(fluid_model_part_name)
            model.AddModelPart(fluid_mesh_model_part)

        super(ALENavierStokesSolverVMSMonolithic, self).__init__(model, navier_stokes_settings)
        # create mesh motion solver
        aux_mesh_solver_settings = KratosMultiphysics.Parameters("""{
            "parallel_type" : "OpenMP"
        }""")
        custom_settings.AddValue("problem_data", aux_mesh_solver_settings)
        custom_settings.AddValue("solver_settings", custom_settings["ale_settings"])
        custom_settings.RemoveValue("ale_settings")

        mesh_motion_solver_settings = custom_settings["solver_settings"]

        if not mesh_motion_solver_settings.Has("echo_level"):
            mesh_motion_solver_settings.AddEmptyValue("echo_level")
            fluid_echo_lvl = navier_stokes_settings["echo_level"].GetInt()
            mesh_motion_solver_settings["echo_level"].SetInt(fluid_echo_lvl)

        if mesh_motion_solver_settings.Has("model_part_name"):
            if not fluid_model_part_name == mesh_motion_solver_settings["model_part_name"].GetString():
                raise Exception('Fluid- and Mesh-Solver have to use the same "model_part_name"!')
        else:
            mesh_motion_solver_settings.AddEmptyValue("model_part_name")
            fluid_model_part_name = navier_stokes_settings["model_part_name"].GetString()
            mesh_motion_solver_settings["model_part_name"].SetString(fluid_model_part_name)

        if mesh_motion_solver_settings.Has("domain_size"):
            fluid_domain_size = navier_stokes_settings["domain_size"].GetInt()
            mesh_motion_domain_size = mesh_motion_solver_settings["domain_size"].GetInt()
            if not fluid_domain_size == mesh_motion_domain_size:
                raise Exception('Fluid- and Mesh-Solver have to use the same "domain_size"!')
        else:
            mesh_motion_solver_settings.AddEmptyValue("domain_size")
            fluid_domain_size = navier_stokes_settings["domain_size"].GetInt()
            mesh_motion_solver_settings["domain_size"].SetInt(fluid_domain_size)

        self.ale_interface_part_names = []
        if mesh_motion_solver_settings.Has("ale_interface_parts"):
            for i in range(mesh_motion_solver_settings["ale_interface_parts"].size()):
                self.ale_interface_part_names.append(
                    mesh_motion_solver_settings["ale_interface_parts"][i].GetString())
            mesh_motion_solver_settings.RemoveValue("ale_interface_parts")

        ## TODO add warning if "compute_reactions" is set to false!

        import python_solvers_wrapper_mesh_motion
        self.ale_solver = python_solvers_wrapper_mesh_motion.CreateSolver(model, custom_settings)

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
        for part_name in self.ale_interface_part_names:
            part_nodes = self.main_model_part.GetSubModelPart(part_name).Nodes
            KratosMultiphysics.VariableUtils().CopyVectorVar(KratosMultiphysics.MESH_VELOCITY,
                                                             KratosMultiphysics.VELOCITY,
                                                             part_nodes)

        self.GetFluidSolver().SolveSolutionStep()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
