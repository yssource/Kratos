import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication
import swimming_DEM_solver
import sympy as sp
import numpy as np
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver
import L2_error_projection_utility as error_projector

class FluidFractionTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        self.project_parameters = project_parameters
        super(FluidFractionTestSolver, self).__init__(model,
                                                      project_parameters,
                                                      field_utility,
                                                      fluid_solver,
                                                      dem_solver,
                                                      variables_manager)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):
        self.ImposeVelocity()
        super(FluidFractionTestSolver, self).SolveFluidSolutionStep()

    def ImposeVelocity(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, 0.0)
            node.Fix(Kratos.VELOCITY_Z)

    def SetUpResultsDatabase(self):
        current_time = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        x = np.array([node.X for node in self.fluid_solver.main_model_part.Nodes])
        y = np.array([node.Y for node in self.fluid_solver.main_model_part.Nodes])
        z = np.array([node.Z for node in self.fluid_solver.main_model_part.Nodes])

        from KratosMultiphysics.SwimmingDEMApplication.custom_body_force import casas_fluid_fraction_solution
        self.casas_solution = casas_fluid_fraction_solution.CasasFluidFractionSolution(self.project_parameters["fluid_parameters"]["processes"]["boundary_conditions_process_list"][1]["Parameters"]["benchmark_parameters"])
        fluid_fraction_list = np.array([self.casas_solution.alpha(current_time, x, y, z) for x, y, z in zip(x, y, z)])
        fluid_fraction_rate_list = np.array([self.casas_solution.dalphat(current_time, x, y, z) for x, y, z in zip(x, y, z)])


        iterator = 0
        for node in self.fluid_solver.main_model_part.Nodes:
            fluid_fraction = fluid_fraction_list[iterator]
            fluid_fraction_rate = fluid_fraction_rate_list[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION, fluid_fraction)
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_RATE, fluid_fraction_rate)
            iterator += 1

    def ConstructL2ErrorProjector(self):
        self.L2_error_projector = error_projector.L2ErrorProjectionUtility(self.fluid_solver.main_model_part)

    def ProjectL2Error(self):
        self.velocity_error_projected, self.pressure_error_projected, self.error_model_part = self.L2_error_projector.ProjectL2()
        return self.velocity_error_projected, self.pressure_error_projected, self.error_model_part

    def SolveDEM(self):
        super(FluidFractionTestSolver, self).SolveDEM()
