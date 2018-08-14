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
        io_settings = data_settings["io_settings"]
        data_array = np.array([])
        cs_tools.ImportArrayFromSolver(from_client, data_name, data_array)

        sdof_solver = self.solvers[self.solver_name]

        value = sum(data_array)

        if "io_options" in io_settings:
            if "swap_sign" in io_settings["io_options"]:
                value *= -1.0

        data_identifier = data_settings["data_identifier"]
        sdof_solver.SetData(data_identifier, value)


    def ExportData(self, data_settings, to_client):
        if not data_settings["data_format"] == "scalar_value":
            raise Exception('SDofIO can only handle scalar values')
        sdof_solver = self.solvers[self.solver_name]

        data_identifier = data_settings["data_identifier"]

        x = sdof_solver.GetData(data_identifier)

        data_settings["scalar_value"] = x

        to_client.ImportData(data_settings, sdof_solver)
