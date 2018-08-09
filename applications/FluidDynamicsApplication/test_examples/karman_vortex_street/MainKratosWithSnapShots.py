from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

from fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time

import numpy as np 
from copy import deepcopy

def GenerateSampleNodeIDs(main_modelpart):
    return [1,2]


class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        self.dimension = 2
        self.fluid_sol_snap_shots = []
        self.step = 0
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            Dt = parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
        
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            self.StoreSolutionSnapShot(model["MainModelPart"])        
            self.step = self.step + 1

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

    def GetSolutionSnapShots(self):
        return self.fluid_sol_snap_shots, self.step

    def StoreSolutionSnapShot(self, main_modelpart):
        current_snap_shot = [0]*self.dimension*len(main_modelpart.Nodes)
        index = 0
        for node in main_modelpart.Nodes:
            current_snap_shot[index*self.dimension + 0] =  node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
            current_snap_shot[index*self.dimension + 1] =  node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,0)
            if(self.dimension == 3):
                current_snap_shot[index*self.dimension + 1] =  node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z,0)
            index = index + 1
            
        self.fluid_sol_snap_shots.append(current_snap_shot)


if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithFlush(model,parameters)
    simulation.Initialize()
    main_modelpart = model["MainModelPart"]
    sampling_nodes_ids = GenerateSampleNodeIDs(main_modelpart)
    snap_shots = []
    dimension = 2
    length_snap_shots = dimension*len(main_modelpart.Nodes)

    for node_id in sampling_nodes_ids:
        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())        
        model = KratosMultiphysics.Model()
        simulation = FluidDynamicsAnalysisWithFlush(model,parameters)
        # Modify the node with the node id 
        # node = main_modelpart.Nodes[node_id]
        # node.X = node.X + 0.1 
        # node.SetSolutionStepValue(VELOCITY_X,0,10.0099)        
        simulation.Run()
        x, steps = simulation.GetSolutionSnapShots()
        snap_shots.append(x)
    
    snap_shots = np.array(snap_shots)
    snap_shots = snap_shots.reshape((len(sampling_nodes_ids)*steps, length_snap_shots))
    snap_shots = np.transpose(snap_shots)
    print(snap_shots.shape)
    