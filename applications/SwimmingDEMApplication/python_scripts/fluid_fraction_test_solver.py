import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication
import swimming_DEM_solver
import sympy as sp
import numpy as np
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver
import L2_error_projection_utility as error_projector

def CartesianToPolar(x, y):
    x0 = 0.0
    y0 = 0.0
    L = 1
    x_rel, y_rel = x-x0, y-y0
    r = np.sqrt(x_rel**2 + y_rel**2)
    theta = x/L
    return r, theta

class FluidFractionTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super(FluidFractionTestSolver, self).__init__(model,
                                                      project_parameters,
                                                      field_utility,
                                                      fluid_solver,
                                                      dem_solver,
                                                      variables_manager)

    def ReturnExactFluidFraction(self, x, y):
        time = self.fluid_solver.main_model_part.ProcessInfo[Kratos.TIME]
        r, theta = CartesianToPolar(x, y)
        alpha0 = 0.7
        alpha_min = 0.5
        assert(alpha0 >= alpha_min and alpha0 <= 1.0 and alpha_min > 0.0)
        delta_alpha = min(alpha0 - alpha_min, 1.0 - alpha0)
        T = 0.1
        omega = 2 * np.pi / T

        return alpha0 + delta_alpha * np.sin(theta + omega * time)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):
        self.ImposeVelocity()
        super(FluidFractionTestSolver, self).SolveFluidSolutionStep()
        self.SetFluidFractionField()

    def SetFluidFractionField(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            fluid_fraction = self.ReturnExactFluidFraction(node.X, node.Y)
            node.SetSolutionStepValue(Kratos.FLUID_FRACTION, fluid_fraction)

    def ImposeVelocity(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, 0.0)
            node.Fix(Kratos.VELOCITY_Z)

    def ConstructL2ErrorProjector(self):
        self.L2_error_projector = error_projector.L2ErrorProjectionUtility(self.fluid_solver.main_model_part)

    def ProjectL2Error(self):
        self.velocity_error_projected, self.pressure_error_projected, self.error_model_part = self.L2_error_projector.ProjectL2()
        return self.velocity_error_projected, self.pressure_error_projected, self.error_model_part

    def SolveDEM(self):
        super(FluidFractionTestSolver, self).SolveDEM()
