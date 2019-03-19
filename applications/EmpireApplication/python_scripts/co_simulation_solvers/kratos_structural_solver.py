from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.ChimeraApplication as kchim

import KratosMultiphysics as km

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from structural_mechanics_analysis import StructuralMechanicsAnalysis

def CreateSolver(cosim_solver_settings, level):
    return KratosStructuralSolver(cosim_solver_settings, level)

class StructuralWithVTKoutput(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super(StructuralWithVTKoutput,self).__init__(model,project_parameters)
        self.parameters = project_parameters

    def Initialize(self):
        super(StructuralWithVTKoutput,self).Initialize()
        main_model_part = self.model["Structure"]
        self.main_model_part = main_model_part
                
       #self.vtkOutput = km.VtkOutput(main_model_part,"nnn",self.parameters["output_configuration"])
        self.step=0

    def OutputSolutionStep(self):
        super(StructuralWithVTKoutput,self).OutputSolutionStep()
        
        #if(self.step%1==0):
        #    self.vtkOutput.PrintOutput()
        #self.step+=1

class KratosStructuralSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return StructuralWithVTKoutput(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__