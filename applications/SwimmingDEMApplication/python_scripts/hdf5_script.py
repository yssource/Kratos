# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.SwimmingDEMApplication as Dem

def Norm(vector):
    return np.sqrt(sum(v**2 for v in vector))

class ErrorProjectionPostProcessTool(object):
    def __init__(self, error_model_part, test_number):

        if not test_number:
            return

        self.parameters = Kratos.Parameters( """
        {
            "file_name": "sp_data.hdf5",
            "target_porosity" : 0.3,
            "probe_height" : 0.2032
        }  """ )

        if test_number == 1: # CTW16
            test_id = "CTW16"
        elif test_number == 2: # CTW10
            test_id = "CTW10"
        else: # Blind test
            test_id = "BlindTest"

        self.parameters.AddEmptyValue("test_id")
        self.parameters["test_id"].SetString(test_id)

        self.error_model_part = error_model_part
        self.time = self.error_model_part.ProcessInfo[Kratos.TIME]
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())

        self.dtype = np.float64
        self.compression_type = 'gzip'
        self.number_of_readings = 0

        self.last_time = - float('inf')

        # This assumes annulus centered at the origin
        #with h5py.File(self.file_path, 'w') as f:
        #     f.attrs['test_id'] = self.parameters["test_id"].GetString()
        #    f.attrs['time'] = self.time
        #     f.attrs['error_velocity'] = self.error_velocity
        #     f.attrs['error_pressure'] = self.error_pressure

    def WriteData(self, velocity_error_projected, pressure_error_projected):
        time = self.time
        name = str(self.number_of_readings + 1)
        with h5py.File(self.file_path, 'w') as f:
            #f.attrs['test_id'] = self.parameters["test_id"].GetString()
            f.create_group(name = name)
            f[name].attrs['test_id'] = "CTWD"
            # f[name].attrs['error_velocity'] = velocity_error_projected
            # f[name].attrs['error_pressure'] = pressure_error_projected
        #column_shape = (len(self.error_model_part.Nodes), )
        #self.ids_array = np.array([node.Id for node in self.error_model_part.Nodes])
        self.v_error = velocity_error_projected
        self.p_error = pressure_error_projected

        if not self.last_time == time:
            with h5py.File(self.file_path, 'r+') as f:
                f.create_dataset(name + '/time', shape = (), dtype = self.dtype)
                f[name + '/time'] = self.time
                f.create_dataset('/velocity_error', shape = (), dtype = self.dtype)
                f['/velocity_error'] = self.v_error
                f.create_dataset('/pressure_error', shape = (), dtype = self.dtype)
                f['/pressure_error'] = self.p_error

        self.last_time = time
        self.number_of_readings += 1