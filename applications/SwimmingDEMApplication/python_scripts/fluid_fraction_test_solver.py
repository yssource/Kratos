import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector
import swimming_DEM_solver
import sympy as sp
import numpy as np
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver

class FluidFractionTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super(FluidFractionTestSolver, self).__init__(model,
                                                      project_parameters,
                                                      field_utility,
                                                      fluid_solver,
                                                      dem_solver,
                                                      variables_manager)

    def ReturnExactFluidFraction(self, x, y):
        #interpolate_process_data = self.project_parameters['processes']['check_interpolated_fluid_fraction'][0]
        #interpolate_process_parameters = interpolate_process_data['Parameters']
        #field_def = [entry.GetString() for entry in interpolate_process_parameters['value']]
        #field = eval('-0.4 * x - 0.4 * y + 1')
        field = 1.0
        return field

    def CannotIgnoreFluidNow(self):
        return self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):
        self.ImposeVelocity()
        #self.SetBodyForceField()
        super(FluidFractionTestSolver, self).SolveFluidSolutionStep()

    def SetFluidFractionField(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            fluid_fraction = self.ReturnExactFluidFraction(node.X, node.Y)
            node.SetSolutionStepValue(Kratos.FLUID_FRACTION, fluid_fraction)

    def ImposeVelocity(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, 0.0)
            node.Fix(Kratos.VELOCITY_Z)

    def SolveDEM(self):
        super(FluidFractionTestSolver, self).SolveDEM()