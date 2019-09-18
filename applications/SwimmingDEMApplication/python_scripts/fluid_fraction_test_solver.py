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
        super(FluidFractionTestSolver, self).__init__(model,
                                                      project_parameters,
                                                      field_utility,
                                                      fluid_solver,
                                                      dem_solver,
                                                      variables_manager)

    def ReturnExactFluidFraction(self, x, y):
        field = eval('-0.4 * x - 0.4 * y + 1')
        #field = 0.5
        return field

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):
        self.ImposeVelocity()
        super(FluidFractionTestSolver, self).SolveFluidSolutionStep()

    def SetFluidFractionField(self):
        #pass
        #field_utilities.PorosityField(0, True).ImposePorosityField(self.fluid_solver.main_model_part)
        for node in self.fluid_solver.main_model_part.Nodes:
            fluid_fraction = self.ReturnExactFluidFraction(node.X, node.Y)
            node.SetSolutionStepValue(Kratos.FLUID_FRACTION, fluid_fraction)
            node.SetSolutionStepValue(Kratos.FLUID_FRACTION_OLD, fluid_fraction)

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
