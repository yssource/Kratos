from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.time_integration_base_scheme import TimeIntegrationBaseScheme
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

        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        # time step
        self.dt = scheme_settings["time_step"]

        # placeholders initial values and predictions
        # initial displacement, velocity and acceleration
        self.u0 = None
        self.v0 = None
        self.a0 = None

        # initial displacement, velocity and acceleration
        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0

        self.force = None

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

    def Predict(self):
        """
        """
        return 2.0 * self.u1 - self.u0

    def _AssembleLHS(self, model):
        """
        """
        return model.m + np.dot(model.b, self.dt) + np.dot(model.k, self.dt ** 2)

    def _AssembleRHS(self, model):
        """
        """
        RHS = self.force - np.dot(model.b, self.v1) - np.dot(model.k, self.u1)  - np.dot(model.k, self.v1 * self.dt)
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