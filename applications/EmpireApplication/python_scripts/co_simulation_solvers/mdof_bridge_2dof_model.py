from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.mdof_solver import MDoFSolver

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
    inertia, stifness and damping are un-coupled
    """
    def __init__(self, cosim_solver_settings, level):

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        '''
        sample json input for the model (system) should be

        "system_parameters":
        {
            "length_of_section" : 0.625,
            "translational_dof" :
            {
                "mass"              : 5.83,
                "target_frequency"  : 2.93,
                "damping_log_decr"  : 0.05
            },
            "rotational_dof" :
            {
                "inertia"           : 0.034,
                "target_frequency"  : 2.64,
                "damping_log_decr"  : 0.107
            }
        }

        maybe initial conditions should be added as well
        '''

        l = parameters["system_parameters"]["length_of_section"]

        # heave (translation)
        # mass over unit length as input
        m = parameters["system_parameters"]["translational_dof"]["mass"]
        m_h = m * l
        f_h = parameters["system_parameters"]["translational_dof"]["target_frequency"]
        k_h = m_h / ((2 * np.pi * f_h) ** 2)
        logd_h = parameters["system_parameters"]["translational_dof"]["damping_ratio"]

        # pitch (rotation)
        # inertia over unit length as input
        I = parameters["system_parameters"]["rotational_dof"]["mass"]
        m_r = I * l
        f_r = parameters["system_parameters"]["rotational_dof"]["target_frequency"]
        k_r = m_r / ((2 * np.pi * f_r) ** 2)
        logd_r = parameters["system_parameters"]["rotational_dof"]["damping_ratio"]

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

    def _GetNodalCoordinates(self, level_height, num_of_levels):
        nodal_coordinates = level_height * np.arange(1,num_of_levels+1)
        return nodal_coordinates

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