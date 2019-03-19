from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CoSimulationBaseConvergenceAccelerator(object):
    def __init__(self, settings, data):
        self.settings = settings
        self.data = data
        self.echo_level = 0

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeCouplingIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.prev_data = self.data.GetNumpyArray()

    def FinalizeCouplingIteration(self):
        pass

    def ComputeUpdate(self):
        current_data = self.data.GetNumpyArray()
        residual = current_data - self.prev_data
        new_data = self.prev_data + self._ComputeUpdate(residual, self.prev_data)
        self.data.ApplyUpdateToData(new_data)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.PrintInfo(("Convergence Accelerator", cs_tools.bcolors.BOLD + self._Name()+cs_tools.bcolors.ENDC))

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        raise Exception('"_Name" has to be implemented in the derived class!')

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def _ComputeUpdate( self, residual, previous_data ):
        raise Exception('"_ComputeUpdate" has to be implemented in the derived class!')
