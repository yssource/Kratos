import KratosMultiphysics
import numpy as np
import sympy as sp

## Import base class file
from custom_body_force.manufactured_solution import ManufacturedSolution

def CreateManufacturedSolution(custom_settings):
    return StationaryFluidFractionSolution(custom_settings)

class StationaryFluidFractionSolution(ManufacturedSolution):
    def __init__(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "velocity"         : 1.0,
                "length"           : 1.0,
                "viscosity"        : 0.1,
                "density"          : 1.0,
                "frequency"        : 1.0,
                "damping"          : 1.0,
                "center_x1"        : 0.0,
                "center_x2"        : 0.0,
                "independent_term" : 0.4,
                "maximum_alpha"    : 1.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.U = settings["velocity"].GetDouble()
        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho
        self.centerx1 = settings["center_x1"].GetDouble()
        self.centerx2 = settings["center_x2"].GetDouble()
        self.independent_term = settings["independent_term"].GetDouble()
        self.maximum_alpha = settings["maximum_alpha"].GetDouble()

    def u1(self, time, x1, x2, x3):
        return  100*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def u2(self, time, x1, x2, x3):
        return 100*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du1dt(self, time, x1, x2, x3):
        return -100*np.pi*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.sin(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) - 100*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du2dt(self, time, x1, x2, x3):
        return -100*np.pi*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.sin(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) - 100*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du11(self, time, x1, x2, x3):
        return 100*self.independent_term*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 100*(-2*self.centerx1 + 2*x1)*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) + 100*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(-2*self.centerx1 + 2*x1 - 2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du12(self, time, x1, x2, x3):
        return 100*self.independent_term*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 100*(-self.centerx1 + x1)**2*(self.centerx1 - x1 + 1)**2*(200*(-2*self.centerx2 + 2*x2)*(-2*self.centerx2 + 2*x2 - 2) + 200*(-self.centerx2 + x2)**2 + 200*(self.centerx2 - x2 + 1)**2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du21(self, time, x1, x2, x3):
        return 100*self.independent_term*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 100*(-self.centerx2 + x2)**2*(self.centerx2 - x2 + 1)**2*(-200*(-2*self.centerx1 + 2*x1)*(-2*self.centerx1 + 2*x1 - 2) - 200*(-self.centerx1 + x1)**2 - 200*(self.centerx1 - x1 + 1)**2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)


    def du22(self, time, x1, x2, x3):
        return 100*self.independent_term*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 100*(-2*self.centerx2 + 2*x2)*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) + 100*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(-2*self.centerx2 + 2*x2 - 2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)


    def du111(self, time, x1, x2, x3):
        return 200*self.independent_term**2*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**3 + 200*self.independent_term*(-2*self.centerx1 + 2*x1)*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 200*self.independent_term*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(-2*self.centerx1 + 2*x1 - 2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 200*(-2*self.centerx1 + 2*x1)*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(-2*self.centerx1 + 2*x1 - 2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) + 200*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) + 200*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du122(self, time, x1, x2, x3):
        return 200*self.independent_term**2*(-self.centerx1 + x1)**2*(100*(-2*self.centerx2 + 2*x2)*(self.centerx2 - x2 + 1)**2 + 100*(-self.centerx2 + x2)**2*(-2*self.centerx2 + 2*x2 - 2))*(self.centerx1 - x1 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**3 + 200*self.independent_term*(-self.centerx1 + x1)**2*(self.centerx1 - x1 + 1)**2*(200*(-2*self.centerx2 + 2*x2)*(-2*self.centerx2 + 2*x2 - 2) + 200*(-self.centerx2 + x2)**2 + 200*(self.centerx2 - x2 + 1)**2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 100*(-self.centerx1 + x1)**2*(self.centerx1 - x1 + 1)**2*(-2400*self.centerx2 + 2400*x2 - 1200)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def du211(self, time, x1, x2, x3):
        return 200*self.independent_term**2*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**3 + 200*self.independent_term*(-self.centerx2 + x2)**2*(self.centerx2 - x2 + 1)**2*(-200*(-2*self.centerx1 + 2*x1)*(-2*self.centerx1 + 2*x1 - 2) - 200*(-self.centerx1 + x1)**2 - 200*(self.centerx1 - x1 + 1)**2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 100*(-self.centerx2 + x2)**2*(2400*self.centerx1 - 2400*x1 + 1200)*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)


    def du222(self, time, x1, x2, x3):
        return 200*self.independent_term**2*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**3 + 200*self.independent_term*(-2*self.centerx2 + 2*x2)*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 200*self.independent_term*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(-2*self.centerx2 + 2*x2 - 2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)**2 + 200*(-2*self.centerx2 + 2*x2)*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(-2*self.centerx2 + 2*x2 - 2)*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) + 200*(-self.centerx2 + x2)**2*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha) + 200*(-100*(-2*self.centerx1 + 2*x1)*(self.centerx1 - x1 + 1)**2 - 100*(-self.centerx1 + x1)**2*(-2*self.centerx1 + 2*x1 - 2))*(self.centerx2 - x2 + 1)**2*np.exp(-time)*np.cos(np.pi*time)/(-self.independent_term*(-self.centerx1 + x1) - self.independent_term*(-self.centerx2 + x2) + self.maximum_alpha)

    def alpha(self, time, x1, x2, x3):
        return -self.independent_term * x1 - self.independent_term * x2 + self.maximum_alpha

    def alpha1(self, time, x1, x2, x3):
        return -self.independent_term

    def alpha2(self, time, x1, x2, x3):
        return -self.independent_term

    def alpha3(self, time, x1, x2, x3):
        return 0.0
