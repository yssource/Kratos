from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
from co_simulation_tools import classprint

def Create(predictor_settings, solvers, cosim_solver_details, level):
    return LinearDerivativeBasedPredictor(predictor_settings, solvers, cosim_solver_details, level)

class LinearDerivativeBasedPredictor(CosimulationBasePredictor):
    def Predict(self):
        data_sizes = [] # saving the sizes of the data to later split them again
        new_data = []
        size_counter = 0
        for data_entry in self.settings["data_list"]:
            data = self._ImportData(data_entry, 1)
            new_data.append(data)
            size_counter += data.size
            data_sizes.append(size_counter)

        combined_new_data = np.concatenate(new_data)

        new_derivative = []
        for data_entry in self.settings["derivative_list"]:
            new_derivative.append(self._ImportData(data_entry, 1))
        combined_new_derivative = np.concatenate(new_derivative)

        #compute linear prediction
        combined_new_data += self.delta_time * combined_new_derivative
        updated_data = np.split(combined_new_data, data_sizes)

        for data_entry, data_update in zip(self.settings["data_list"], updated_data):
            self._ExportData(data_entry, data_update)

        if self.echo_level > 3:
            classprint(self.lvl, self._Name(), "Computed prediction")

    def _Name(self):
        return self.__class__.__name__

