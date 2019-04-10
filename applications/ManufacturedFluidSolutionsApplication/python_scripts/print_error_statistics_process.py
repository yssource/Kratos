# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ManufacturedFluidSolutionsApplication as MS

# Other imports
import h5py
import time
import numpy as np

from KratosMultiphysics.ManufacturedFluidSolutionsApplication.manufactured_solution_base_process import ManufacturedSolutionBaseProcess as ManufacturedProcess

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PrintErrorStatisticsProcess(Model, settings["Parameters"])

## All the python manufactured processes should be derived from a base class, which is derived from "Process"
class PrintErrorStatisticsProcess(ManufacturedProcess):
    '''
    A wrapper for all the processes related to the manufactured solution,
    such as the body force, the boundary conditions o the error computation processes
    '''
    def __init__(self, model, settings ):

        default_settings = KM.Parameters("""
            {
                "model_part_name"  : "model_part",
                "file_name"        : "output_file",
                "manufactured_id"  : 0,
                "label"            : "description"    
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.f = h5py.File(self.settings["file_name"].GetString() + ".hdf5", 'a') # 'a' means append mode

        dset_name = "manufactured_{:03d}".format(self.settings["manufactured_id"].GetInt())
        if dset_name in self.f:
            self.dset = self.f[dset_name]
            self.dataset_is_initialized = True
        else:
            case_dtype = np.dtype([("label", h5py.special_dtype(vlen=str)),
                                ("num_nodes", np.uint32),
                                ("num_elems", np.uint32),
                                ("time_step", np.float),
                                ("computational_time", np.float),
                                ("rel_error", np.float)])
            self.dset = self.f.create_dataset(dset_name, (0,), maxshape=(None,), chunks=True, dtype=case_dtype)
            self.dataset_is_initialized = False

        self.start_time = time.time()


    def ExecuteFinalize(self):
        self._WriteBodyForceAttributes()
        self._WriteAverageRelativeError()


    def _WriteBodyForceAttributes(self):
        manufactured_parameters = self.manufactured_solution.GetParameters()
        for attr, param in manufactured_parameters.items():
            if param.IsBool():
                value = param.GetBool()
            if param.IsInt():
                value = param.GetInt()
            if param.IsDouble():
                value = param.GetDouble()
            if self.dataset_is_initialized:
                if self.dset.attrs[attr] != value:
                    self.dset.attrs["Warning"] = "There are several vortex definitions in this dataset"
            self.dset.attrs[attr] = value


    def _WriteAverageRelativeError(self):
        elapsed_time = time.time() - self.start_time
        rel_err = self.manufactured_process.ComputeMean(MS.VELOCITY_RELATIVE_ERROR)
        case_data = (self.settings["label"].GetString(),
                     self.model_part.NumberOfNodes(),
                     self.model_part.NumberOfElements(),
                     self.model_part.ProcessInfo[KM.DELTA_TIME],
                     elapsed_time,
                     rel_err)

        case_idx = self.dset.len()
        self.dset.resize((case_idx+1,))
        self.dset[case_idx] = case_data
