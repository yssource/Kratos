from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import *
import KratosMultiphysics.StructuralMechanicsApplication

from structural_mechanics_analysis import StructuralMechanicsAnalysis

import math

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

class StructuralMechanicsAnalysisWithCentrifugalForces(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super(StructuralMechanicsAnalysisWithCentrifugalForces,self).__init__(model,project_parameters)
        self.frequency = 1.0
        
    def Initialize(self):
        super(StructuralMechanicsAnalysisWithCentrifugalForces,self).Initialize()


    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        self.main_model_part = model["Structure"]
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            Dt = parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
            
            self.ApplyCentrifugalBCs(self.main_model_part)
            self.InitializeSolutionStep()           
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()            
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def ApplyCentrifugalBCs(self, model_part ):

        for node in model_part.Nodes:

            radius_vector = Vector(3)
            radius_vector[0] = node.X
            radius_vector[1] = node.Y
            radius_vector[2] = 0
            radius = radius_vector[0]*radius_vector[0] + radius_vector[1]*radius_vector[1] + radius_vector[2]*radius_vector[2]

            normalized_radius_vector = Vector(3)
            if(radius!=0):
                normalized_radius_vector[0] = radius_vector[0]/radius
                normalized_radius_vector[1] = radius_vector[1]/radius
                normalized_radius_vector[2] = radius_vector[2]/radius
            else:
                normalized_radius_vector[0] = 0.0
                normalized_radius_vector[1] = 0.0
                normalized_radius_vector[2] = 0.0 

            omega = 5.0 * math.sin(self.frequency*self.time)

            centrigual_acc = Vector(3)
            centrigual_acc = omega*omega * radius * normalized_radius_vector

            node.SetSolutionStepValue(VOLUME_ACCELERATION,0,centrigual_acc)

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysisWithCentrifugalForces(model,parameters)
    simulation.Run()
