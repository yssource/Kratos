import KratosMultiphysics
import numpy as np
import sympy as sp

## Import base class file
from custom_body_force.manufactured_solution import ManufacturedSolution

def CreateManufacturedSolution(custom_settings):
    return StationaryFluidFractionVortex(custom_settings)

class StationaryFluidFractionVortex(ManufacturedSolution):
    def __init__(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "velocity"    : 1.0,
                "length"      : 1.0,
                "viscosity"   : 0.1,
                "density"     : 1.0,
                "frequency"   : 1.0,
                "damping"     : 1.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.U = settings["velocity"].GetDouble()
        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho

    def f(self, x):
        return 100 * x**2 * (self.L - x)**2

    def g(self, t):
        return sp.cos(np.pi * t) * sp.exp(-t)

    def df(self, x):
        return 2 * 100 * x * (self.L - x)**2 - 2 * 100 * x**2 * (self.L - x)

    def dg(self, t):
        return -np.pi * np.sin(np.pi * t) * np.e**(-t) -np.cos(np.pi * t) * np.e**(-t)

    def ddf(self, x):
        return 2 * 100 * (self.L**2 - 6*self.L*x + 6*x**2)

    def dddf(self, x):
        return 12 * 100 * (2*x - self.L)

    def u1(self, t, x1, x2, x3):
        return self.f(x1) * self.df(x2) * self.g(t)

    def u2(self, t, x1, x2, x3):
        return -self.df(x1) * self.f(x2) * self.g(t)

    def du1dt(self, t, x1, x2, x3):
        return self.f(x1) * self.df(x2) * self.dg(t)

    def du2dt(self, t, x1, x2, x3):
        return -self.df(x1) * self.f(x2) * self.dg(t)

    def du11(self, t, x1, x2, x3):
        return self.df(x1) * self.df(x2) * self.g(t)

    def du12(self, t, x1, x2, x3):
        return self.f(x1) * self.ddf(x2) *  self.g(t)

    def du21(self, t, x1, x2, x3):
        return -self.ddf(x1) * self.f(x2) * self.g(t)

    def du22(self, t, x1, x2, x3):
        return -self.df(x1) * self.df(x2) * self.g(t)

    def du111(self, t, x1, x2, x3):
        return self.ddf(x1) * self.df(x2) * self.g(t)

    def du122(self, t, x1, x2, x3):
        return self.f(x1) * self.dddf(x2) * self.g(t)

    def du211(self, t, x1, x2, x3):
        return -self.dddf(x1) * self.f(x2) * self.g(t)

    def du222(self, t, x1, x2, x3):
        return -self.df(x1) * self.ddf(x2) * self.g(t)

    def alpha(self, t, x1, x2, x3):
        self.fluid_fraction = -0.4 * x1 -0.4 * x2 + 1
        #self.fluid_fraction = 1.0
        return self.fluid_fraction


