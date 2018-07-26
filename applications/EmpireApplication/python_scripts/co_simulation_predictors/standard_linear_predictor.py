from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from co_simulation_base_predictor import CosimulationBasePredictor

import numpy as np
import co_simulation_ios.co_simulation_io_factory as io_factory

from co_simulation_tools import csprint, blue


def Create(predictor_settings, solvers, cosim_solver_details, level):
    return StandardLinearPredictor(predictor_settings, solvers, cosim_solver_details, level)


class StandardLinearPredictor(CosimulationBasePredictor):
    # @param beta factor for weighting last and current value of the predicted values. Can be set in interval: [0, 1.0]
    def __init__(self, settings, solvers, cosim_solver_details, level):
        super().__init__(settings, solvers, cosim_solver_details, level)    #Warum wird hier bei den anderen Klassen das Argument übergeben? Das ist nicht notwendig!
        self.settings = settings
        self.solvers = solvers
        self.cosim_solver_details = cosim_solver_details
        self.io = io_factory.CreateIO(settings, solvers, "None", cosim_solver_details, level)
        self.lvl = level

    def Predict(self):
        data_sizes = [] # saving the sizes of the data to later split them again
        size_counter = 0
        new_data = []
        old_data = []
        for data_entry in self.settings["data_list"]:
            data = self._ImportData(data_entry, 0)    #besser?  data = super(AveragedTractionPredictor, self)._ImportData(data_entry, 0)
            new_data.append(data)
            size_counter += data.size
            data_sizes.append(size_counter)
        self.combined_new_data = np.concatenate(new_data)

        for data_entry in self.settings["old_data_list"]:
            old_data.append(self._ImportData(data_entry, 1))
        combined_old_data = np.concatenate(old_data)

        #compute prediction
        self.combined_new_data = 2 * self.combined_new_data - combined_old_data
        updated_data = np.split(self.combined_new_data, data_sizes)

        for data_entry, data_update in zip(self.settings["data_list"], updated_data):
            self._ExportData(data_entry, data_update)

        csprint(self.lvl, blue("<< Compute prediction with standard_linear_predictor >>"))

    def _Name(self):
        return self.__class__.__name__

