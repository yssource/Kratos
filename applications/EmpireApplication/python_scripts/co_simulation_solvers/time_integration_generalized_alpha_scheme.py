from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.time_integration_base_scheme import TimeIntegrationBaseScheme
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateScheme(scheme_settings):
    return TimeIntegrationGeneralizedAlphaScheme(scheme_settings)

class TimeIntegrationGeneralizedAlphaScheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "type" : "generalized_alpha",
                "time_step" : 0.01,
                "settings": {
                    "p_inf" : 0.15
                }
            }

        RecursivelyValidateAndAssignDefaults(default_settings, scheme_settings)

        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        # time step
        self.dt = scheme_settings["time_step"]
        pInf = scheme_settings["settings"]["p_inf"]
        self.alphaM = (2.0 * pInf - 1.0) / (pInf + 1.0)
        self.alphaF = pInf / (pInf + 1.0)
        self.beta = 0.25 * (1 - self.alphaM + self.alphaF)**2
        self.gamma = 0.5 - self.alphaM + self.alphaF

        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.dt**2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alphaF

        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)

        #coefficients for damping
        self.a1b = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a2b = (1.0 - self.alphaF) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alphaF) * (0.5 * self.gamma / self.beta - 1.0) * self.dt

        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF

        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt

        # coefficients for acceleration update
        self.a1a = self.a1v / (self.dt * self.gamma)
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)

        # placeholders initial values and predictions
        # initial displacement, velocity and acceleration
        self.u0 = None
        self.v0 = None
        self.a0 = None

        # initial displacement, velocity and acceleration
        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0

		# force from a previous time step (initial force)
        self.f0 = None
        self.f1 = None

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

		# force from a previous time step (initial force)
        self.f0 = np.dot(model.m,self.a0) + np.dot(model.b,self.v0) + np.dot(model.k,self.u0)
        self.f1 = np.dot(model.m,self.a1) + np.dot(model.b,self.v1) + np.dot(model.k,self.u1)

    def Predict(self):
        """
        """
        return 2.0 * self.u1 - self.u0

    def _AssembleLHS(self, model):
        """
        """
        return self.a1h * model.m + self.a2h * model.b + self.a3h * model.k

    def _AssembleRHS(self, model):
        """
        """
        # PMT -> F = ... -> self.f1 is actually for this step a new external force which is read in
        # so use SetSolutionStep() in case there is an external force
        # should be AssembleRHS(self, model, f1)
        #F = (1.0 - self.alphaF) * f1 + self.alphaF * self.f0

        f = (1.0 - self.alphaF) * self.force + self.alphaF * self.f0
        RHS = np.dot(model.m,(self.a1m * self.u0 + self.a2m * self.v0 + self.a3m * self.a0))
        RHS += np.dot(model.b,(self.a1b * self.u0 + self.a2b * self.v0 + self.a3b * self.a0))
        RHS += np.dot(self.a1k * model.k, self.u0) + f

        # and here an update of self.f1
        self.f1 = self.force
        return RHS

    def UpdateDerivedValues(self):
        """
        """
        self.v1 = self.a1v * (self.u1 - self.u0) + self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + self.a2a * self.v0 + self.a3a * self.a0