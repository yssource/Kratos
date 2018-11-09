# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from compute_drag_process import ComputeDragProcess

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeDragProcessForChimera(model, settings["Parameters"])

class ComputeDragProcessForChimera(ComputeDragProcess):
    def __init__(self, model, params ):
        super(ComputeDragProcessForChimera,self).__init__(model,params)
        self.previoustime = 0
    
    def _GetFileHeader(self):
        header  = '# Drag for model part ' + self.params["model_part_name"].GetString() + '\n'
        header += '# Time Fx Fy Fz\n'
        return header
    
    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if time != self.previoustime:
            super(ComputeDragProcessForChimera,self).ExecuteFinalizeSolutionStep()
            self.previoustime = time
    
    def _GetCorrespondingDragForce(self):
        return KratosCFD.DragUtilities().CalculateBodyFittedDrag(self.model_part)