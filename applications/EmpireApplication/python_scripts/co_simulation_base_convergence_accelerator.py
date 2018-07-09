from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import co_simulation_io_factory
import numpy as np

class CoSimulationBaseConvergenceAccelerator(object):
    def __init__(self, settings, solvers, cosim_solver_details, level):
        self.settings = settings
        self.solvers = solvers
        self.cosim_solver_details = cosim_solver_details
        self.lvl = level
        self.echo_level = 0
        if "echo_level" in self.settings:
            self.echo_level = self.settings["echo_level"]
        self.io = co_simulation_io_factory.CreateIO(settings, solvers, "None", cosim_solver_details, level)

    def AdvanceInTime(self):
        # Saving the previous data for the computation of the residual
        previous_data = [] # discard previous data fields
        self.data_sizes = []
        size_counter = 0
        for data_entry in self.settings["data_list"]:
            prev_data = self.__ImportData(data_entry)
            previous_data.append(prev_data)
            size_counter += prev_data.size
            self.data_sizes.append(size_counter)

        self.combined_prev_data = np.concatenate(previous_data)

    def ComputeUpdate(self):
        new_data = []
        for data_entry in self.settings["data_list"]:
            new_data.append(self.__ImportData(data_entry))

        combined_residuals = np.concatenate(new_data) - self.combined_prev_data

        combined_data_update = self._ComputeUpdate(combined_residuals, self.combined_prev_data)

        data_updates = np.split(combined_data_update, self.data_sizes)

        for data_entry, data_update in zip(self.settings["data_list"], data_updates):
            self.__ExportData(data_entry, data_update)

    def ImportData(self, data_name, from_client):
        pass
    def ImportMesh(self, mesh_name, from_client):
        pass

    def ExportData(self, data_name, to_client):
        pass
    def ExportMesh(self, mesh_name, to_client):
        pass

    def MakeDataAvailable(self, data_name, to_client):
        pass
    def MakeMeshAvailable(self, mesh_name, to_client):
        pass

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
