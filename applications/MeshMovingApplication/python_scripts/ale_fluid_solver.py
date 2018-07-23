from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication

# Other imports
from python_solver import PythonSolver
import python_solvers_wrapper_fluid


def CreateSolver(model, custom_settings):
    '''This function creates the requested solver
    If no "ale_settings" are specified a regular fluid-solver is created
    '''
    if custom_settings.Has("ale_settings"):
        ALEFluidSolver(model, custom_settings)
    else:
        KratosMultiphysics.PrintInfo("ALEFluidSolver", "No ale settings found, creating a pure fluid solver")
        return python_solvers_wrapper_fluid.CreateSolver(model, custom_settings)


class ALEFluidSolver(PythonSolver):
    def __init__(self, model, custom_settings):
        super(ALEFluidSolver, self).__init__(model, custom_settings["solver_settings"])

        ## Creating the fluid solver
        self.fluid_solver = self._CreateFluidSolver(custom_settings)

        ## Creating the mesh-motion solver
        KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")
        mesh_motion_settings = custom_settings.Clone()
        # delete the (fluid) solver settings
        mesh_motion_settings.RemoveValue("solver_settings")
        # add the ale solver settings as solver settings
        mesh_motion_settings.AddValue("solver_settings", custom_settings["ale_settings"])

        import python_solvers_wrapper_mesh_motion
        self.mesh_motion_solver = python_solvers_wrapper_mesh_motion.CreateSolver(model, mesh_motion_settings)

        # Getting the min_buffer_size from both solvers
        # and assigning it to the fluid_solver, bcs this one handles the model_part
        self.fluid_solver.min_buffer_size = max(self.fluid_solver.GetMinimumBufferSize(),
                                                self.mesh_motion_solver.GetMinimumBufferSize())

        self.is_printing_rank = self.fluid_solver._IsPrintingRank()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Construction finished")

    def AddVariables(self):
        self.fluid_solver.AddVariables()
        self.mesh_motion_solver.AddVariables()
        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Variables Added")

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.mesh_motion_solver.AddDofs()
        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "DOFs Added")

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.mesh_motion_solver.Initialize()
        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("::[ALEFluidSolver]::", "Finished initialization")

    def ImportModelPart(self):
        self.fluid_solver.ImportModelPart() # only ONE solver imports the ModelPart

    def PrepareModelPart(self):
        # Doing it ONLY for the fluid solver
        return self.fluid_solver.PrepareModelPart()

    def AdvanceInTime(self, current_time):
        # Doing it ONLY for the fluid solver
        return self.fluid_solver.AdvanceInTime(current_time)

    def Finalize(self):
        self.fluid_solver.Finalize()
        self.mesh_motion_solver.Finalize()

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.mesh_motion_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.mesh_motion_solver.Predict()

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.mesh_motion_solver.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        self.mesh_motion_solver.SolveSolutionStep()

    def Check(self):
        self.fluid_solver.Check()
        self.mesh_motion_solver.Check()

    def Clear(self):
        self.fluid_solver.Clear()
        self.mesh_motion_solver.Clear()


    def GetFluidSolver(self):
        return self.fluid_solver

    def GetMeshMotionSolver(self):
        return self.mesh_motion_solver

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()


    def _CreateFluidSolver(self, custom_settings):
        '''This function creates the fluid solver.
        It can be overridden to create different fluid solvers
        '''
        fluid_settings = custom_settings.Clone()
        # remove the ale_settings so we can use the fluid_solver_wrapper constructor
        fluid_settings["solver_settings"].RemoveValue("ale_settings")
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolver(model, fluid_settings)
