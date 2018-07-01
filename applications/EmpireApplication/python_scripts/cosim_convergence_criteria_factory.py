import numpy as np
from numpy import linalg as la
import math


def Create(settings):
    pass


class CoSimulationConvergenceCriteria(object):
    def __init__(self, settings):
        pass

    def AdvanceInTime(self):
        pass

    def IsConverged(self, data_field):
        return la.norm(data_field)/math.sqrt(data_field.size) < self.abs_tol