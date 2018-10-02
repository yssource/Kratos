from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.mdof_solver import MDoFSolver
from co_simulation_tools import ValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateSolver(cosim_solver_settings, level):
    return MDoFBridge2DoFModel(cosim_solver_settings, level)

class MDoFBridge2DoFModel(MDoFSolver):
    """
    MDoF model - 2DoF for a bridge section

    ATTENTION:
    For this model it is assumed that
    inertia, stiffness and damping are decoupled
    """
    def __init__(self, cosim_solver_settings, level):

        input_file_name = cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        default_settings = json.loads("""{
                "system_parameters":{
                    "length_of_section" : 0.625,
                    "help"              : "1st value - translational dof, 2nd value - rotational dof",
                    "mass_per_length"   : [5.83, 0.034],
                    "target_frequency"  : [2.93, 2.64],
                    "damping_log_decr"  : [0.05, 0.107]
                },
                "initial_conditions":{
                    "displacement"      : [0.1, 0.025],
                    "velocity"          : [0.0, 0.0],
                    "acceleration"      : [0.0, 0.0],
                    "external_load"     : [0.0, 0.0]
                },
                "time_integration_parameters":{},
                "solver_parameters":{},
                "output_parameters":{}
            }""")

        ValidateAndAssignDefaults(default_settings, parameters)

        l = parameters["system_parameters"]["length_of_section"]

        # heave (translation)
        # mass over unit length as input

        m = parameters["system_parameters"]["mass_per_length"][0]
        m_h = m * l
        f_h = parameters["system_parameters"]["target_frequency"][0]
        k_h = m_h * ((2 * np.pi * f_h) ** 2)
        logd_h = parameters["system_parameters"]["damping_log_decr"][0]

        # pitch (rotation)
        # inertia over unit length as input
        I = parameters["system_parameters"]["mass_per_length"][1]
        m_r = I * l
        f_r = parameters["system_parameters"]["target_frequency"][1]
        k_r = m_r * ((2 * np.pi * f_r) ** 2)
        logd_r = parameters["system_parameters"]["damping_log_decr"][1]

        m = self._CalculateMass(m_h, m_r)
        k = self._CalculateStiffness(k_h, k_r)
        b = self._CalculateDamping(m_h, m_r,
                                   k_h, k_r,
                                   logd_h, logd_r)

        # creating a model dictionary to pass to base class constructor
        model = {}
        model.update({'M': m})
        model.update({'K': k})
        model.update({'B': b})

        initial_conditions = self._SetupInitialValues(parameters['initial_conditions'])
        model.update({'initial_conditions': initial_conditions})

        super(MDoFBridge2DoFModel, self).__init__(model, cosim_solver_settings, level)

    def _CalculateMass(self, m_h, m_r):
        """
        Calculate mass m
        """
        return np.array([[m_h, 0],
                         [0, m_r]])

    def _CalculateStiffness(self, k_h, k_r):
        """
        Calculate mass k
        """
        return np.array([[k_h, 0],
                         [0, k_r]])

    def _CalculateDamping(self, m_h, m_r, k_h, k_r, logd_h, logd_r):
        """
        Calculate damping b
        """
        xi_h = logd_h / (2 * np.pi)
        xi_r = logd_r / (2 * np.pi)

        return np.array([[2 * xi_h * np.sqrt(m_h * k_h), 0],
                         [0, 2 * xi_r * np.sqrt(m_r * k_r)]])

    def _GetIOName(self):
        return "mdof_bridge_2dof_model"

    def _Name(self):
        return self.__class__.__name__

    # PMT: to be implemented
    def _DofList(self):
        '''
        A DoF list saying which DoF entry
        what kind of deformation it represents
        In this case probably:
        ["DeltaY","ThethaZ"]
        '''
        pass

    def _SetupInitialValues(self, initial_values):
        '''
        From a list generate numpy array for compatibility
        '''
        for key, value in initial_values.items():

            value = np.array([num_val for num_val in value], dtype='float64')

            initial_values[key] = value

        return initial_values