from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

import KratosMultiphysics.ChimeraApplication as kchim

try:
    import KratosMultiphysics.MeshMovingApplication
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "not imported")

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidDynamicsAnalysisWithVTK(FluidDynamicsAnalysis):
    def __init__(self,model,project_parameters):
        super(FluidDynamicsAnalysisWithVTK,self).__init__(model,project_parameters)

    def Initialize(self):
        super(FluidDynamicsAnalysisWithVTK,self).Initialize()
        main_model_part = self.model["FluidModelPart"]

        fluid = main_model_part.GetSubModelPart("Parts_Fluid")
        self.VtkOut =kchim.VtkOutput(fluid,"nnn",self.project_parameters["output_configuration"])
        self.step = 0


    def OutputSolutionStep(self):
        super(FluidDynamicsAnalysisWithVTK,self).OutputSolutionStep()
        if(self.step%100==0):
            self.VtkOut.PrintOutput()
        self.step+=1


def CreateSolver(cosim_solver_settings, level):
    return KratosFluidSolver(cosim_solver_settings, level)

class KratosFluidSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysisWithVTK(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__

