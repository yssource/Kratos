from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
from co_simulation_tools import classprint

def Create(predictor_settings, solvers, cosim_solver_details, level):
    return StandardLinearPredictor(predictor_settings, solvers, cosim_solver_details, level)

class StandardLinearPredictor(CosimulationBasePredictor):
    def __init__(self, settings, solvers, cosim_solver_details, level):
        super(StandardLinearPredictor, self).__init__(settings, solvers, cosim_solver_details, level)
        # TODO add comment why we do this
        num_data = len(self.settings["data_list"])
        self.data_arrays_t1 = [np.array([]) for e in range(num_data)]
        self.data_arrays_t2 = [np.array([]) for e in range(num_data)]

        # TODO check buffer size!

    def Predict(self):
        for i, data_entry in enumerate(self.settings["data_list"]):
            solver = self.solvers[data_entry["solver"]]
            data_name = data_entry["data_name"]
            cs_tools.ImportArrayFromSolver(solver, data_name, self.data_arrays_t1[i], 1) # should be 0/1 ?
            cs_tools.ImportArrayFromSolver(solver, data_name, self.data_arrays_t2[i], 2) # should be 1/2 ?

            self.data_arrays_t1[i] = 2*self.data_arrays_t1[i] - data_arrays_t2[i]

        self._UpdateData(self.data_arrays_t1)

    def _Name(self):
        return self.__class__.__name__

