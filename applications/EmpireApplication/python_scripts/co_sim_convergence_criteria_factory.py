from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import io_factory

import numpy as np
from numpy import linalg as la
import math

from co_simulation_tools import csprint, bold, green, red


def CreateConvergenceCriteria(settings, solvers, cosim_solver_details, level):
    return CoSimulationConvergenceCriteria(settings, solvers, cosim_solver_details, level)

class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings, solvers, cosim_solver_details, level):
        self.settings = settings
        self.io = io_factory.CreateIO(settings, solvers, "None", cosim_solver_details, level)
        self.solvers = solvers
        self.echo_level = 0
        if "echo_level" in self.settings:
            self.echo_level = self.settings["echo_level"]
        self.lvl = level
        self.abs_tolerances = []
        self.rel_tolerances = []
        for data_entry in self.settings["data_list"]:
            self.abs_tolerances.append(data_entry["abs_tolerance"])
            self.rel_tolerances.append(data_entry["rel_tolerance"])
        data_size = len(self.settings["data_list"])

    def AdvanceInTime(self):
        self.old_data = [] # discard old data fields
        for data_entry in self.settings["data_list"]:
            self.old_data.append(self.__ImportData(data_entry))

    def IsConverged(self):
        convergence_list = []
        idx = 0
        for data_entry in self.settings["data_list"]:
            new_data = self.__ImportData(data_entry)
            residual = new_data - self.old_data[idx]
            abs_norm = la.norm(residual) / math.sqrt(residual.size)
            norm_new_data = la.norm(new_data)
            if norm_new_data < 1e-15:
                norm_new_data = 1.0 # to avoid division by zero
            rel_norm = la.norm(residual) / norm_new_data
            convergence_list.append(abs_norm < self.abs_tolerances[idx] or rel_norm < self.rel_tolerances[idx])
            if self.echo_level > 0:
                info_msg  = 'Convergence for "'+bold(data_entry["data_name"])+'": '
                if convergence_list[idx]:
                    info_msg += green("ACHIEVED")
                else:
                    info_msg += red("NOT ACHIEVED")
                csprint(self.lvl, info_msg)
            if self.echo_level > 1:
                info_msg  = bold("abs_norm")+" = " + str(abs_norm) + " | "
                info_msg += bold("abs_tol")+" = " + str(self.abs_tolerances[idx])
                info_msg += " || "+bold("rel_norm")+" = " + str(rel_norm) + " | "
                info_msg += bold("rel_tol") +" = " + str(self.rel_tolerances[idx])
                csprint(self.lvl, info_msg)
            idx += 1

        return min(convergence_list) # return false if any of them did not converge!

    def __ImportData(self, data_entry):
        data_name = data_entry["data_name"]
        from_solver = data_entry["from_solver"]
        return self.io.ImportData(data_name, self.solvers[from_solver])