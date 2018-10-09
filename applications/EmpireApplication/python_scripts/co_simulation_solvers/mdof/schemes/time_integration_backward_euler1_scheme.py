from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from time_integration_base_scheme import TimeIntegrationBaseScheme
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateScheme(scheme_settings):
    return TimeIntegrationBackwardEuler1Scheme(scheme_settings)

# PMT to be checked, seems to be written acceleration based, might need correction

class TimeIntegrationBackwardEuler1Scheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "type" : "backward_euler1",
                "time_step" : 0.01,
            }

        RecursivelyValidateAndAssignDefaults(default_settings, scheme_settings)

        super(TimeIntegrationBackwardEuler1Scheme, self).__init__(scheme_settings)

    def _AssembleLHS(self, model):
        """
        """
        return model.m

    def _AssembleRHS(self, model):
        """
        """
        self.f0 = self.force

        RHS = self.f1 * self.dt**2
        RHS -= np.dot(- 2 * model.m + model.b * self.dt, self.u1)
        RHS -= np.dot(model.m - model.b * self.dt + model.k * self.dt**2, self.u2)

        return RHS

    def UpdateDerivedValues(self):
        """
        """
        self.v0 = (self.u0 - self.u1) / self.dt
        self.a0 = (self.u0 - 2 * self.u1 + self.u2) / self.dt**2