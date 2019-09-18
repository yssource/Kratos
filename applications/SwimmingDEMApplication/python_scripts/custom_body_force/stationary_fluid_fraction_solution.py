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
        # self.X = sp.symbols("X")
        # self.T  = sp.symbols("T")
        # self.fx = 100 * self.X**2 * (self.L - self.X)**2
        # self.gt = sp.cos(np.pi * self.T) * sp.exp(-self.T)
        # self.dfx = sp.diff(self.fx, self.X)
        # self.dgt = sp.diff(self.gt, self.T)
        # self.ddfx = sp.diff(self.fx, self.X, 2)
        # self.dddfx = sp.diff(self.fx, self.X, 3)

    def f(self, x):
        return 100 * x**2 * (self.L - x)**2

    def g(self, t):
        return np.cos(np.pi * t) * np.e**(-t)

    def df(self, x):
        return 2 * 100 * x * (self.L - x)**2 - 2 * 100 * x**2 * (self.L - x)

    def dg(self, t):
        return -np.pi * np.sin(np.pi * t) * np.e**(-t) - np.cos(np.pi * t) * np.e**(-t)

    def ddf(self, x):
        return 2 * 100 * (self.L**2 - 6 * self.L * x + 6 * x**2)

    def dddf(self, x):
        return 12 * 100 * (2 * x - self.L)

    def u1(self, t, x1, x2, x3):
        return self.f(x1) * self.df(x2) * self.g(t) / self.alpha(t, x1, x2, x3)

    def u2(self, t, x1, x2, x3):
        return -self.df(x1) * self.f(x2) * self.g(t) / self.alpha(t, x1, x2, x3)

    def du1dt(self, t, x1, x2, x3):
        return self.f(x1) * self.df(x2) * self.dg(t) / self.alpha(t, x1, x2, x3)

    def du2dt(self, t, x1, x2, x3):
        return -self.df(x1) * self.f(x2) * self.dg(t) / self.alpha(t, x1, x2, x3)

    def du11(self, t, x1, x2, x3):
        return (self.df(x1) * self.df(x2) * self.g(t) * self.alpha(t, x1, x2, x3) - self.f(x1) * self.df(x2) * self.g(t) * self.dalpha1(x1, x2)) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du12(self, t, x1, x2, x3):
        return (self.f(x1) * self.ddf(x2) *  self.g(t) * self.alpha(t, x1, x2, x3) - self.f(x1) * self.df(x2) *  self.g(t) * self.dalpha2(x1, x2)) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du21(self, t, x1, x2, x3):
        return (-self.ddf(x1) * self.f(x2) * self.g(t) * self.alpha(t, x1, x2, x3) - (-self.df(x1) * self.f(x2) * self.g(t) * self.dalpha1(x1, x2))) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du22(self, t, x1, x2, x3):
        return (-self.df(x1) * self.df(x2) * self.g(t) * self.alpha(t, x1, x2, x3) - (-self.df(x1) * self.f(x2) * self.g(t) * self.dalpha2(x1, x2))) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du111(self, t, x1, x2, x3):
        return ((self.ddf(x1) * self.df(x2) * self.g(t) * self.alpha(t, x1, x2, x3) + self.df(x1) * self.df(x2) * self.g(t) * self.dalpha1(x1, x2) - (self.df(x1) * self.df(x2) * self.g(t) * self.dalpha1(x1, x2) + self.f(x1) * self.df(x2) * self.g(t) * self.ddalpha1(x1, x2))) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) - (self.df(x1) * self.df(x2) * self.g(t) * self.alpha(t, x1, x2, x3) - self.f(x1) * self.df(x2) * self.g(t) * self.dalpha1(x1, x2)) * (self.dalpha1(x1, x2) * self.alpha(t, x1, x2, x3) + self.alpha(t, x1, x2, x3) * self.dalpha1(x1, x2))) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du122(self, t, x1, x2, x3):
        return ((((self.f(x1) * self.dddf(x2) * self.g(t) * self.alpha(t, x1, x2, x3) + self.f(x1) * self.ddf(x2) *  self.g(t) * self.dalpha2(x1, x2)) - (self.f(x1) * self.ddf(x2) *  self.g(t) * self.dalpha2(x1, x2) + self.f(x1) * self.df(x2) *  self.g(t) * self.ddalpha2(x1, x2))) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3)) - (self.f(x1) * self.ddf(x2) *  self.g(t) * self.alpha(t, x1, x2, x3) - self.f(x1) * self.df(x2) *  self.g(t) * self.dalpha2(x1, x2)) * (self.dalpha2(x1, x2) * self.alpha(t, x1, x2, x3) + self.alpha(t, x1, x2, x3) * self.dalpha2(x1, x2))) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du211(self, t, x1, x2, x3):
        return (((((-self.dddf(x1) * self.f(x2) * self.g(t) * self.alpha(t, x1, x2, x3) + (-self.ddf(x1) * self.f(x2) * self.g(t) * self.dalpha1(x1, x2))) - ((-self.ddf(x1) * self.f(x2) * self.g(t) * self.dalpha1(x1, x2) + (-self.df(x1) * self.f(x2) * self.g(t) * self.ddalpha1(x1, x2)))))) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3)) - (-self.ddf(x1) * self.f(x2) * self.g(t) * self.alpha(t, x1, x2, x3) - (-self.df(x1) * self.f(x2) * self.g(t) * self.dalpha1(x1, x2))) * (self.dalpha1(x1, x2) * self.alpha(t, x1, x2, x3) + self.alpha(t, x1, x2, x3) * self.dalpha1(x1, x2))) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def du222(self, t, x1, x2, x3):
        return ((((-self.df(x1) * self.ddf(x2) * self.g(t) * self.alpha(t, x1, x2, x3) + (-self.df(x1) * self.df(x2) * self.g(t) * self.dalpha2(x1, x2))) - (-self.df(x1) * self.df(x2) * self.g(t) * self.dalpha2(x1, x2) + (-self.df(x1) * self.f(x2) * self.g(t) * self.ddalpha2(x1, x2)))) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3)) - ((-self.df(x1) * self.df(x2) * self.g(t) * self.alpha(t, x1, x2, x3) - (-self.df(x1) * self.f(x2) * self.g(t) * self.dalpha2(x1, x2)))) * (self.dalpha2(x1, x2) * self.alpha(t, x1, x2, x3) + self.alpha(t, x1, x2, x3) * self.dalpha2(x1, x2))) / (self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3) * self.alpha(t, x1, x2, x3))

    def alpha(self, t, x1, x2, x3):
        self.fluid_fraction = -0.4 * x1 - 0.4 * x2 + 1
        return self.fluid_fraction

    def dalpha1(self, x1, x2):
        self.dfluid_fraction = -0.4
        return self.dfluid_fraction

    def dalpha2(self, x1, x2):
        self.dfluid_fraction = -0.4
        return self.dfluid_fraction

    def ddalpha1(self, x1, x2):
        self.ddfluid_fraction = 0.0
        return self.ddfluid_fraction

    def ddalpha2(self, x1, x2):
        self.ddfluid_fraction = 0.0
        return self.ddfluid_fraction

