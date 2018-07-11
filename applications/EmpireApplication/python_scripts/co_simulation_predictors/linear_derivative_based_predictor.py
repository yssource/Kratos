from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7


def Create(predictor_settings, solvers, cosim_solver_details, level):
    return LinearDerivativeBasedPredictor(predictor_settings, solvers, cosim_solver_details, level)

class LinearDerivativeBasedPredictor(object):
    def __init__(self, settings, solvers, cosim_solver_details, level):
        self.settings = settings
        self.solvers = solvers
        self.cosim_solver_details = cosim_solver_details
        self.lvl = level
        self.echo_level = 0
        if "echo_level" in self.settings:
            self.echo_level = self.settings["echo_level"]

    def Predict(self):
        print("Hier könnte ihre Werbung stehen")