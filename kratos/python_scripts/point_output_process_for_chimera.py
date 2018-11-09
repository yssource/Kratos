from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from point_output_process import PointOutputProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PointOutputProcessForChimera(Model, settings["Parameters"])

class PointOutputProcessForChimera(PointOutputProcess):
    def __init__(self, model, params):
        super(PointOutputProcessForChimera,self).__init__(model,params)
        self.previoustime = 0
    
    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if time != self.previoustime:
            super(PointOutputProcessForChimera,self).ExecuteFinalizeSolutionStep()
            self.previoustime = time
        
