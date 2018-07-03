from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import io_factory

import numpy as np
from numpy import linalg as la
import math

from co_simulation_tools import csprint, bold, green, red


def CreateConvergenceCriteria(settings, solvers, level):
    return CoSimulationConvergenceCriteria(settings, solvers, level)

class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings, solvers, level):
        self.settings = settings
        self.io = io_factory.CreateIO(settings)
        self.solvers = solvers
        self.echo_level = 0
        if "echo_level" in self.settings:
            self.echo_level = self.settings["echo_level"]
        self.lvl = level
        self.tolerances = []
        for data_entry in self.settings["data_list"]:
            self.tolerances.append(data_entry["tolerance"])
        data_size = len(self.settings["data_list"])

    def InitializeSolutionStep(self):
        self.old_data = [] # discard old data fields
        for data_entry in self.settings["data_list"]:
            self.old_data.append(self.__ImportData(data_entry))

    def IsConverged(self):
        convergence_list = []
        idx = 0
        for data_entry in self.settings["data_list"]:
            residual = self.__ImportData(data_entry) - self.old_data[idx]
            norm = la.norm(residual) / math.sqrt(residual.size)
            convergence_list.append(norm < self.tolerances[idx])
            if self.echo_level > 0:
                info_msg  = 'Convergence for "'+bold(data_entry["data_name"])+'": '
                if convergence_list[idx]:
                    info_msg += green("ACHIEVED")
                else:
                    info_msg += red("NOT ACHIEVED")
                info_msg += " : norm = " + str(norm) + " | tol = " + str(self.tolerances[idx])
                csprint(self.lvl, info_msg)
            idx += 1

        return min(convergence_list) # return false if any of them did not converge!

    def __ImportData(self, data_entry):
        data_name = data_entry["data_name"]
        geometry_name = data_entry["geometry_name"]
        from_solver = data_entry["from_solver"]
        return self.io.ImportData(data_name, geometry_name, self.solvers[from_solver])