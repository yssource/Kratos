from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import co_simulation_ios.co_simulation_io_factory as io_factory
import numpy as np
from co_simulation_tools import classprint, bold

class CosimulationBasePredictor(object):
    def __init__(self, settings, solvers, cosim_solver_details, level):
        self.settings = settings
        self.solvers = solvers
        self.cosim_solver_details = cosim_solver_details
        self.lvl = level
        self.echo_level = 0
        self.io = io_factory.CreateIO(settings, solvers, "None", cosim_solver_details, level)

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def SetDeltaTime(self, delta_time):
        self.delta_time = delta_time

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        classprint(self.lvl, "Predictor", bold(self._Name()))

    def Check(self):
        print("The predictors do not yet implement Check!")

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        raise Exception('"_Name" has to be implemented in the derived class!')


    #ATTENTION: Problem with private functions in classes:
    # _ private method -> when calling it in subclass _functionname
    # __very private method -> when calling it in subclass classname__function name might cause problems

    def _ImportData(self, data_entry, buffer_index = 0):
        data_name = data_entry["data_name"]
        from_solver = data_entry["from_solver"]
        data_definition = self.solvers[from_solver].GetDataDefinition(data_name)
        old_buffer_index = 0
        if "buffer_index" in data_definition:
            old_buffer_index = data_definition["buffer_index"]
        data_definition["buffer_index"] = buffer_index
        data = self.io.ImportData(data_name, self.solvers[from_solver])
        data_definition["buffer_index"] = old_buffer_index
        return data

    def _ExportData(self, data_entry, data_update):
        data_name = data_entry["data_name"]
        from_solver = data_entry["from_solver"] # This is "from_solver", bcs we give back the updated solution
        return self.io.ExportData(data_name, self.solvers[from_solver], data_update)

