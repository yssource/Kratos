from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import co_simulation_ios.co_simulation_io_factory as io_factory
import numpy as np
from co_simulation_tools import classprint, bold

class CoSimulationBaseConvergenceAccelerator(object):
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

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        previous_data = [] # discard previous data fields
        self.data_sizes = [] # saving the sizes of the data to later split them again
        size_counter = 0
        for data_entry in self.settings["data_list"]:
            prev_data = self.__ImportData(data_entry)
            previous_data.append(prev_data)
            size_counter += prev_data.size
            self.data_sizes.append(size_counter)

        self.combined_prev_data = np.concatenate(previous_data)

    def FinalizeNonLinearIteration(self):
        pass

    def ComputeUpdate(self):
        new_data = []
        for data_entry in self.settings["data_list"]:
            new_data.append(self.__ImportData(data_entry))

        combined_residuals = np.concatenate(new_data) - self.combined_prev_data

        combined_new_data = self.combined_prev_data + self._ComputeUpdate(combined_residuals, self.combined_prev_data)

        updated_data = np.split(combined_new_data, self.data_sizes)

        for data_entry, data_update in zip(self.settings["data_list"], updated_data):
            self.__ExportData(data_entry, data_update)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        classprint(self.lvl, "Convergence Accelerator", bold(self._Name()))

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        raise Exception('"_Name" has to be implemented in the derived class!')

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def __ImportData(self, data_entry):
        data_name = data_entry["data_name"]
        from_solver = data_entry["from_solver"]
        return self.io.ImportData(data_name, self.solvers[from_solver])

    def __ExportData(self, data_entry, data_update):
        data_name = data_entry["data_name"]
        from_solver = data_entry["from_solver"] # This is "from_solver", bcs we give back the updated solution
        return self.io.ExportData(data_name, self.solvers[from_solver], data_update)

    def _ComputeUpdate( self, residual, previous_data ):
        raise Exception('"_ComputeUpdate" has to be implemented in the derived class!')
