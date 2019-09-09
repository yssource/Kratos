from KratosMultiphysics import Model, Parameters, Logger
import swimming_DEM_procedures as SDP
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication
import os
import sys
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
sys.path.insert(0, dir_path)
from swimming_DEM_analysis import SwimmingDEMAnalysis
from swimming_DEM_analysis import Say
import L2_error_projection_utility as error_projector
import fluid_fraction_test_solver as sdem_solver

class FluidFractionTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, iteration, varying_parameters = Parameters("{}")):

        from KratosMultiphysics.SwimmingDEMApplication import hdf5_script
        self.projector_post_process = hdf5_script.ErrorProjectionPostProcessTool(iteration)
        super(FluidFractionTestAnalysis, self).__init__(model, varying_parameters)
        self.project_parameters = varying_parameters
        self.iteration = iteration

    def Initialize(self):
        super(FluidFractionTestAnalysis, self).Initialize()
        self._GetSolver().SetFluidFractionField()
        # self._GetSolver().ImposePressure()
        self._GetSolver().ConstructL2ErrorProjector()

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def _CreateSolver(self):
        import fluid_fraction_test_solver as sdem_solver
        return sdem_solver.FluidFractionTestSolver(self.model,
                                                   self.project_parameters,
                                                   self.GetFieldUtility(),
                                                   self._GetFluidAnalysis()._GetSolver(),
                                                   self._GetDEMAnalysis()._GetSolver(),
                                                   self.vars_man)

    def FinalizeSolutionStep(self):
        super(FluidFractionTestAnalysis, self).FinalizeSolutionStep()
        self.velocity_error_projected, self.pressure_error_projected, self.error_model_part = self._GetSolver().ProjectL2Error()
        from KratosMultiphysics.SwimmingDEMApplication import hdf5_script
        self.projector_post_process.WriteData(self.error_model_part, self.velocity_error_projected, self.pressure_error_projected)

    def TransferBodyForceFromDisperseToFluid(self):
        pass

if __name__ == "__main__":
    # Setting parameters

    with open('ProjectParameters.json','r') as parameter_file:
        parameters = Parameters(parameter_file.read())

    # Create Model
    model = Model()

    # To avoid too many prints
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    test = FluidFractionTestAnalysis(model, parameters)
    test.Run()