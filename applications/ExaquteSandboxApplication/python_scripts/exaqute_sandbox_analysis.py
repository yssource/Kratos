from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from sys import argv

import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage
import KratosMultiphysics.ExaquteSandboxApplication

class ExaquteSandboxAnalysis(AnalysisStage):
    '''Analysis Stage only meant to run cases of Exaqute application.'''

    def __init__(self,model,parameters):
        super(ExaquteSandboxAnalysis,self).__init__(model,parameters)

    def ModifyInitialProperties(self):
        model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        for node in self.model.GetModelPart(model_part_name).Nodes:
            coord_x = node.X
            coord_y = node.Y
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing)

    def _CreateSolver(self):
        import KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_stationary_solver
        return KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])

    def _GetSimulationName(self):
        return "Exaqute Sandbox Analysis"

if __name__ == '__main__':
    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "problem_settings/project_parameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ExaquteSandboxAnalysis(model,parameters)
    simulation.Run()
