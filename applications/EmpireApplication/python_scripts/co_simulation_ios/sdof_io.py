from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_base_io import CoSimulationBaseIO

# Other imports
import numpy as np
import co_simulation_tools as cs_tools

def Create(io_settings, solvers, solver_name, cosim_solver_details, level):
    return SDofIO(io_settings, solvers, solver_name, cosim_solver_details, level)

class SDofIO(CoSimulationBaseIO):

    def ImportData(self, data_settings, from_client):
        data_name = data_settings["data_name"]
        data_array = np.array([])
        cs_tools.ImportArrayFromSolver(from_client, data_name, data_array)

        sdof_solver = self.solvers[self.solver_name]

        # Do sth with the data_array, sum etc ...?

        sdof_solver.SetRHS(data_array) # this has to be somehow available, maybe through public interface?


    def ExportData(self, data_settings, to_client):
        sdof_solver = self.solvers[self.solver_name]

        dx = sdof_solver.GetDx() # this has to be somehow available, maybe through public interface?

        data_settings["single_value"] = dx

        to_client.ImportData(data_settings, sdof_solver)
