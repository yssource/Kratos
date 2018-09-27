from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.mdof_solver import MDoFSolver

# Other imports
import numpy as np
import json
import os

def CreateSolver(cosim_solver_settings, level):
    return MDoFSDoFModel(cosim_solver_settings, level)

class MDoFSDoFModel(MDoFSolver):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, cosim_solver_settings, level):

        input_file_name = cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        '''
        sample json input for the model (system) should be

        "system_parameters":
        {
            "mass"              : 1.0,
            "target_frequency"  : 1.5,
            "damping_ratio"     : 0.05,
            "level_height"      : 3.5,
            "number_of_levels"  : 1
        }

        maybe initial conditions should be added as well
        '''

        m = parameters["system_parameters"]["mass"]
        target_freq = parameters["system_parameters"]["target_frequency"]
        zeta = parameters["system_parameters"]["damping_ratio"]

        k = self._CalculateStiffness(m, target_freq)
        b = self._CalculateDamping(m, k, zeta)

        level_height = parameters["system_parameters"]["level_height"]
        num_of_levels = parameters["system_parameters"]["number_of_levels"]

        height_coordinates = self._GetNodalCoordinates(level_height, num_of_levels)

        nodal_coordinates = {"x0": np.zeros(len(height_coordinates)),
                            "y0": height_coordinates,
                            "x": None,
                            "y": None}

        # creating a model dictionary to pass to base class constructor
        model = {}
        model.update({'M': np.array([[m]])})
        model.update({'K': np.array([[k]])})
        model.update({'B': np.array([[b]])})
        # check if this is needed
        model.update({'nodal_coordinates': nodal_coordinates})

        model.update({'initial_values':
                        {'displacement' : np.array([parameters["initial_values"]["displacement"]]),
                         'velocity'     : np.array([parameters["initial_values"]["velocity"]]),
                         'acceleration' : np.array([parameters["initial_values"]["acceleration"]]),
                         'external_load' : np.array([parameters["initial_values"]["external_load"]])
                         }})

        super(MDoFSDoFModel, self).__init__(model, cosim_solver_settings, level)

    def _CalculateStiffness(self, m, target_freq):
        """
        Calculate stiffness k
        """
        return m * (target_freq * 2 * np.pi)**2

    def _CalculateDamping(self, m, k, zeta):
        """
        Calculate damping b
        """
        return zeta * 2.0 * np.sqrt(m * k)

    def _GetNodalCoordinates(self, level_height, num_of_levels):
        nodal_coordinates = level_height * np.arange(1,num_of_levels+1)
        return nodal_coordinates

    def _GetIOName(self):
        return "mdof_sdof_model"

    def _Name(self):
        return self.__class__.__name__

    # PMT: to be implemented
    def _DofList(self):
        '''
        A DoF list saying which DoF entry
        what kind of deformation it represents
        In this case probably:
        ["DeltaX"]
        '''
        pass