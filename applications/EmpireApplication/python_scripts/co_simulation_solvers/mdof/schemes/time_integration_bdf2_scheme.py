from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from time_integration_base_scheme import TimeIntegrationBaseScheme
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateScheme(scheme_settings):
    return TimeIntegrationBDF2Scheme(scheme_settings)

# PMT to be checked, seems to be written acceleration based, might need correction

class TimeIntegrationBDF2Scheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "type" : "bdf2",
                "time_step" : 0.01,
                "settings" : {
                    "nr_of_dofs"    : 1,
                    "buffer_size"   : 4
                }
            }

        RecursivelyValidateAndAssignDefaults(default_settings, scheme_settings)

        # base scheme settings
        super(TimeIntegrationBDF2Scheme, self).__init__(scheme_settings)

    def Initialize(self, model):
        """
        """
        # initial displacement, velocity and acceleration
        self.u0 = model.u0
        self.v0 = model.v0
        self.a0 = model.a0
        # initial force
        # PMT is this needed/correct like this?
        self.force = model.f0

        # initial displacement, velocity and acceleration
        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0

    def _AssembleLHS(self, model):
        """
        """
        return model.m

    def _AssembleRHS(self, model):
        """
        """
        RHS = self.force - np.dot(model.b, self.v0) - np.dot(model.k, self.u0)
        # PMT might be a problem, as it seems to be solving for accelerations
        #self.an0 =  np.linalg.solve(LHS, RHS)
        # should be displacement based
        # also check indexing u0 vs u1 which is new-old

        # and here an update of self.f1
        self.f1 = self.force
        return RHS

    def UpdateDerivedValues(self):
        """
        """
        self.u0 = self.u1 + self.dt * self.v1
        self.v0 = self.v1 + self.dt * self.a1