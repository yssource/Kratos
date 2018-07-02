from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import io_factory

import numpy as np
from numpy import linalg as la
import math


def CreateConvergenceCriteria(settings, solvers):
    return CoSimulationConvergenceCriteria(settings, solvers)

class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings, solvers):
        self.settings = settings
        self.io = io_factory.CreateIO(settings)
        self.solvers = solvers
        self.echo_level = 0

    def IsConverged(self):
        return False
        # # import the data fields
        # for data in self.settings[data_list]:
        #     data_name = data["data_name"]
        #     geometry_name = data["geometry_name"]
        #     from_solver = data["from_solver"]
        #     data_field = self.io.ImportData(data_name, geometry_name, self.solvers[from_solver])

        #     # compute norms and stuff ...

        # return la.norm(data_field)/math.sqrt(data_field.size) < self.abs_tol