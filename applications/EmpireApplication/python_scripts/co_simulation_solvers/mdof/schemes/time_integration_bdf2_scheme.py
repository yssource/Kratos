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

    def _AssembleLHS(self, model):
        """
        """
        pass

    def _AssembleRHS(self, model):
        """
        """
        pass

    def UpdateDerivedValues(self):
        """
        """
        pass