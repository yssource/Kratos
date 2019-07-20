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
                "label"            : "description",
                "variables_list"   : [],
                "write_reynolds"   : true,
                "write_strouhal"   : true
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.f = h5py.File(self.settings["file_name"].GetString() + ".hdf5", 'a') # 'a' means append mode

        self.variables = [KM.KratosGlobals.GetVariable(v) for v in self.settings["variables_list"].GetStringArray()]


    def ExecuteBeforeSolutionLoop(self):
        KM.FindNodalHNonHistoricalProcess(self.model_part).Execute()
        self.dset = self._GetManufacturedDataset()
        self.start_time = time.time()


    def ExecuteFinalizeSolutionStep(self):
        self.manufactured_process.ComputeNodalCFL()


    def ExecuteFinalize(self):
        self._WriteAverageError()


    def _FillAdditionalAttributes(self):
        if self.settings["write_reynolds"].GetBool():
            reynolds = self.manufactured_solution.Reynolds()
            self.manufactured_solution.GetParameters().AddEmptyValue("reynolds").SetDouble(reynolds)
        if self.settings["write_strouhal"].GetBool():
            strouhal = self.manufactured_solution.Strouhal()
            self.manufactured_solution.GetParameters().AddEmptyValue("strouhal").SetDouble(strouhal)
        

    def _WriteAttributes(self, dset):
        manufactured_parameters = self.manufactured_solution.GetParameters()
        for attr, param in manufactured_parameters.items():
            if param.IsBool():
                value = param.GetBool()
            if param.IsInt():
                value = param.GetInt()
            if param.IsDouble():
                value = param.GetDouble()
            dset.attrs[attr] = value


    def _CheckAttributes(self, attributes):
        match = None

        manufactured_parameters = self.manufactured_solution.GetParameters()
        for key, param in manufactured_parameters.items():
            match = False
            if param.IsBool():
                value = param.GetBool()
            if param.IsInt():
                value = param.GetInt()
            if param.IsDouble():
                value = param.GetDouble()

            if attributes.get(key) is not None:
                if attributes[key] == value:
                    match = True
            if not match:
                break

        return match


    def _GetManufacturedDataset(self):
        self._FillAdditionalAttributes()
        for name, data in self.f.items():
            if self._CheckAttributes(data.attrs):
                return data

        dset_name = "manufactured_{:03d}".format(self.f.items().__len__())
        header = self._GetHeaderDtype()
        dset = self.f.create_dataset(dset_name, (0,), maxshape=(None,), chunks=True, dtype=header)
        self._WriteAttributes(dset)
        return dset


    def _GetHeaderDtype(self):
        header = [("label", h5py.special_dtype(vlen=str)),
                  ("num_nodes", np.uint32),
                  ("num_elems", np.uint32),
                  ("time_step", np.float),
                  ("computational_time", np.float)]

        for variable in self.variables:
            header.append((variable.Name(), np.float))
        return np.dtype(header)


    def _WriteAverageError(self):
        elapsed_time = time.time() - self.start_time
        case_data = [self.settings["label"].GetString(),
                     self.model_part.NumberOfNodes(),
                     self.model_part.NumberOfElements(),
                     self.model_part.ProcessInfo[KM.DELTA_TIME],
                     elapsed_time]

        for variable in self.variables:
            value = self.manufactured_process.ComputeMean(variable)
            case_data.append(value)

        case_idx = self.dset.len()
        self.dset.resize((case_idx+1,))
        self.dset[case_idx] = tuple(case_data)
