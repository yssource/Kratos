from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ChimeraApplication
import KratosMultiphysics.ChimeraApplication as kchim
try:
    import KratosMultiphysics.MeshMovingApplication
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "not imported")

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from fluid_chimera_analysis import FluidChimeraAnalysis

def CreateSolver(cosim_solver_settings, level):
    return KratosChimeraSolver(cosim_solver_settings, level)

class ChimeraWithVTKoutput(FluidChimeraAnalysis):

    def __init__(self,model,project_parameters):
        super(ChimeraWithVTKoutput,self).__init__(model,project_parameters)
        
        #output_post  = project_parameters.Has("output_configuration")

    def Initialize(self):
        super(ChimeraWithVTKoutput,self).Initialize()
        main_model_part = self.model["FluidModelPart"]
        self.main_model_part = main_model_part
        background = main_model_part.GetSubModelPart("GENERIC_background")
        self.patch= main_model_part.GetSubModelPart("GENERIC_patch")
        #self.structure = main_model_part.GetSubModelPart("GENERIC_structure")
        
        self.vtkOutput_background = kchim.VtkOutput(background,"nnn",self.parameters["output_configuration"])
        self.vtkOutput_patch = kchim.VtkOutput(self.patch,"nnn",self.parameters["output_configuration"])
        #self.vtkOutput_structure = kchim.VtkOutput(self.structure,"nnn",self.parameters["output_configuration"])
        self.step=0
    def RunSolutionLoop(self):
        """
        this function never called !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        So need to add patch movemnt inside Co simulation analysis , call Initialise solution step here and translate or rotate patch there
        """
        
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            Dt = self.parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            
    def OutputSolutionStep(self):
        super(ChimeraWithVTKoutput,self).OutputSolutionStep()
        
        if(self.step%1==0):
            self.vtkOutput_background.PrintOutput()
            self.vtkOutput_patch.PrintOutput()
        self.step+=1


class KratosChimeraSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return ChimeraWithVTKoutput(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()
    
    def SolveSolutionStep(self):
        self._GetAnalysisStage().FinalizeSolutionStep()
        #self._GetAnalysisStage().ChimeraProcess.ExecuteFinalizeSolutionStep()
        self._GetAnalysisStage().InitializeSolutionStep()
        #self._GetAnalysisStage().ChimeraProcess.ExecuteInitializeSolutionStep()
        self._GetAnalysisStage()._GetSolver().SolveSolutionStep()
        #self._GetAnalysisStage().FinalizeSolutionStep()

    def _Name(self):
        return self.__class__.__name__
