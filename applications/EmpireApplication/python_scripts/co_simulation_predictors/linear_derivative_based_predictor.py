from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7


import co_simulation_ios.co_simulation_io_factory as io_factory

import numpy as np
from numpy import linalg as la

from co_simulation_tools import csprint, bold, green, red, magenta

def Create(predictor_settings, solvers, cosim_solver_details, level):
    return LinearDerivativeBasedPredictor(predictor_settings, solvers, cosim_solver_details, level)

class LinearDerivativeBasedPredictor(object):
    def __init__(self, settings, solvers, cosim_solver_details, level):
        self.settings = settings
        self.solvers = solvers
        self.cosim_solver_details = cosim_solver_details
        self.io = io_factory.CreateIO(settings, solvers, "None", cosim_solver_details, level)
        self.lvl = level
        self.echo_level = 0
        if "echo_level" in self.settings:
            self.echo_level = self.settings["echo_level"]

    def SetDeltaTime(self, delta_time):
        self.delta_time = delta_time

    def Predict(self):

        data_sizes = [] # saving the sizes of the data to later split them again
        new_data = []
        size_counter = 0
        for data_entry in self.settings["data_list"]:
            data = self.__ImportData(data_entry, 1)
            new_data.append(data)
            size_counter += data.size
            data_sizes.append(size_counter)

        combined_new_data = np.concatenate(new_data)

        new_derivative = []
        for data_entry in self.settings["derivative_list"]:
            new_derivative.append(self.__ImportData(data_entry, 1))
        combined_new_derivative = np.concatenate(new_derivative)

        #compute linear prediction
        combined_new_data += self.delta_time * combined_new_derivative
        updated_data = np.split(combined_new_data, data_sizes)

        for data_entry, data_update in zip(self.settings["data_list"], updated_data):
            self.__ExportData(data_entry, data_update)

        csprint(self.lvl, magenta("<< Compute prediction >>"))


    def __ImportData(self, data_entry, buffer_index = 0):
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


    def __ExportData(self, data_entry, data_update):
        data_name = data_entry["data_name"]
        from_solver = data_entry["from_solver"] # This is "from_solver", bcs we give back the updated solution
        return self.io.ExportData(data_name, self.solvers[from_solver], data_update)