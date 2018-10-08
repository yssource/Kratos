from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import json
import os

class TimeIntegrationBaseScheme(object):
    """
    """
    def __init__(self, scheme_settings):

        # provide and validate default setttings
        default_settings = None

        # time step
        self.dt = None

        # placeholders initial values and predictions
        # initial displacement, velocity and acceleration

        # u0 -> u(current-)0
        self.u0 = None
        self.v0 = None
        self.a0 = None

        # initial displacement, velocity and acceleration
        # u1 -> u(current-)1
        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0

        # u2 -> u(current-)2
        self.u2 = self.u1
        self.v2 = self.v1
        self.a2 = self.a1

		# force from a previous time step (initial force)
        self.f0 = None
        self.f1 = None
        self.f2 = None

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

        self.u2 = self.u1
        self.v2 = self.v1
        self.a2 = self.a1

        self.f0 = self.force

    def Predict(self):
        """
        """
        return 2.0 * self.u0 - self.u1

    def _AssembleLHS(self, model):
        """
        """
        pass

    def _AssembleRHS(self, model):
        """
        """
        pass

    def Solve(self, model):
        # sys of eq reads: LHS * u1 = RHS
        LHS = self._AssembleLHS(model)
        RHS = self._AssembleRHS(model)
        self.u0 = np.linalg.solve(LHS, RHS)

    def UpdateDerivedValues(self):
        """
        """
        pass

    # PMT check indices
    # u0 -> current
    # u1 -> old? rather new? being solved for
    def GetDisplacement(self):
        """
        """
        return self.u0

    def GetVelocity(self):
        """
        """
        return self.v0

    def GetAcceleration(self):
        """
        """
        return self.a0

    def GetLoad(self):
        """
        """
        return self.f0

    def GetPreviousDisplacement(self):
        """
        """
        return self.u1

    def GetPreviousVelocity(self):
        """
        """
        return self.v1

    def GetPreviousAcceleration(self):
        """
        """
        return self.a1

    def GetPreviousLoad(self):
        """
        """
        return self.f1

    def AdvanceScheme(self):

        self.u2 = self.u1
        self.v2 = self.v1
        self.a2 = self.a1
        self.f2 = self.f1

        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0
        self.f1 = self.f0