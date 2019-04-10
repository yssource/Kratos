from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication.python_solvers_wrapper_mesh_motion as mesh_mothion_solvers

# Importing the base class
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

class PotentialFlowFSIAnalysis(PotentialFlowAnalysis):
    '''Main script for potential flow simulations with basic fluid-structure-interaction.'''

    def __init__(self,model,project_parameters):
        super(PotentialFlowFSIAnalysis, self).__init__(model,project_parameters)
        self.iteration_count = 0
        self.disp = 0
        self.prev_disp = 0


#    def Run(self):
 #
  #      super(PotentialFlowFSIAnalysis, self).Run()

    def _CreateSolver(self):
        self.mesh_solver = mesh_mothion_solvers.CreateSolver(self.model, self.project_parameters)
        return super(PotentialFlowFSIAnalysis, self)._CreateSolver()

    def StructuralModelResponse(self, lift):
        k = 50 # spring stiffness
        return lift/k


    def KeepAdvancingSolutionLoop(self):
        f = 0.01 # relative change for convergence check
        max_iteration_count = 4
        self.iteration_count = self.iteration_count+1
        if self.disp==0.0:
            return self.iteration_count<max_iteration_count
        else:
            rel_change = abs(1-self.prev_disp/self.disp)
            return self.iteration_count<max_iteration_count and rel_change>f

    def MoveMesh(self):
        self.body_model_part = self.model["MainModelPart.Body2D_BodySurface"]
        for node in self.body_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0,self.disp)


    def RunSolutionLoop(self):
        while self.KeepAdvancingSolutionLoop():
            super(PotentialFlowFSIAnalysis, self).RunSolutionLoop()
            self.prev_disp = self.disp
            self.disp = self.StructuralModelResponse(KratosMultiphysics.TOTAL_LIFT)
            self.MoveMesh()
