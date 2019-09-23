import KratosMultiphysics
import numpy as np
import sympy as sp
from KratosMultiphysics.SwimmingDEMApplication import field_utilities

## Import base class file
from custom_body_force.manufactured_solution import ManufacturedSolution

def CreateManufacturedSolution(custom_settings):
    return TimeDependantFluidFractionVortex(custom_settings)

class TimeDependantFluidFractionVortex(ManufacturedSolution):
    def __init__(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "velocity"    : 1.0,
                "length"      : 1.0,
                "viscosity"   : 0.1,
                "density"     : 1.0,
                "frequency"   : 1.0,
                "damping"     : 1.0,
                "alpha0"      : 0.7,
                "alpha_min"   : 0.5,
                "period"      : 0.01
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.U = settings["velocity"].GetDouble()
        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho
        self.alpha0 = settings["alpha0"].GetDouble()
        period = settings["period"].GetDouble()
        alpha_min = settings["alpha_min"].GetDouble()
        self.delta_alpha = min(self.alpha0 - alpha_min, 1.0 - self.alpha0)
        self.omega = 2 * np.pi / period

    def u1(self, time, x1, x2, x3):
        return (x1*(-self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time))*np.exp(-time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L)))

    def u2(self, time, x1, x2, x3):
        return (-x2*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time))*np.exp(-time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L)))

    def du1dt(self, time, x1, x2, x3):
        return (x1*(self.delta_alpha*self.omega*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) - 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time))*np.cos(self.omega*time + x2/self.L) + (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(self.delta_alpha*self.omega**2*np.exp(time)*np.sin(self.omega*time + x2/self.L) - 20000*np.pi*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.sin(np.pi*time) - 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time)))*np.exp(-time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))**2)

    def du2dt(self, time, x1, x2, x3):
        return x2*(self.delta_alpha*self.omega*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time))*np.cos(self.omega*time + x2/self.L) + (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(self.delta_alpha*self.omega**2*np.exp(time)*np.sin(self.omega*time + x2/self.L) + 20000*np.pi*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.sin(np.pi*time) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time)))*np.exp(-time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))**2

    def du11(self, time, x1, x2, x3):
        return ((-self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time) + 20000*x1*x2*(x1 - 1)*(3*x1 - 1)*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time))*np.exp(-time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L)))

    def du111(self, time, x1, x2, x3):
        return (20000*x2*(x2 - 1)*(2*x2 - 1)*(5*x1*(x1 - 1) + x1*(3*x1 - 1) + (x1 - 1)**2 + (x1 - 1)*(3*x1 - 1))*np.exp(-time)*np.cos(np.pi*time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L)))

    def du12(self, time, x1, x2, x3):
        return (x1*(self.delta_alpha*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) - 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time))*np.cos(self.omega*time + x2/self.L) + (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(20000*self.L*x1*(x1 - 1)**2*(2*x2*(x2 - 1) + x2*(2*x2 - 1) + (x2 - 1)*(2*x2 - 1))*np.cos(np.pi*time) + self.delta_alpha*self.omega*np.exp(time)*np.sin(self.omega*time + x2/self.L)))*np.exp(-time)/(self.L*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))**2))

    def du122(self, time, x1, x2, x3):
        return (-x1*(2*self.delta_alpha*(self.delta_alpha*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) - 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time))*np.cos(self.omega*time + x2/self.L) + (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(20000*self.L*x1*(x1 - 1)**2*(2*x2*(x2 - 1) + x2*(2*x2 - 1) + (x2 - 1)*(2*x2 - 1))*np.cos(np.pi*time) + self.delta_alpha*self.omega*np.exp(time)*np.sin(self.omega*time + x2/self.L)))*np.cos(self.omega*time + x2/self.L) + (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(self.delta_alpha*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) - 20000*x1*x2*(x1 - 1)**2*(x2 - 1)*(2*x2 - 1)*np.cos(np.pi*time))*np.sin(self.omega*time + x2/self.L) - (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(120000*self.L**2*x1*(x1 - 1)**2*(2*x2 - 1)*np.cos(np.pi*time) + self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L))))*np.exp(-time)/(self.L**2*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))**3))

    def du22(self, time, x1, x2, x3):
        return ((-self.L*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time)) + self.delta_alpha*x2*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time))*np.cos(self.omega*time + x2/self.L) - x2*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(20000*self.L*x1*(x1 - 1)*(2*x1 - 1)*(x2 - 1)*(3*x2 - 1)*np.cos(np.pi*time) - self.delta_alpha*self.omega*np.exp(time)*np.sin(self.omega*time + x2/self.L)))*np.exp(-time)/(self.L*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))**2))

    def du222(self, time, x1, x2, x3):
        return ((2*self.delta_alpha*(self.L*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time)) - self.delta_alpha*x2*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time))*np.cos(self.omega*time + x2/self.L) + x2*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(20000*self.L*x1*(x1 - 1)*(2*x1 - 1)*(x2 - 1)*(3*x2 - 1)*np.cos(np.pi*time) - self.delta_alpha*self.omega*np.exp(time)*np.sin(self.omega*time + x2/self.L)))*np.cos(self.omega*time + x2/self.L) - (self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(2*self.L*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(20000*self.L*x1*(x1 - 1)*(2*x1 - 1)*(x2 - 1)*(3*x2 - 1)*np.cos(np.pi*time) - self.delta_alpha*self.omega*np.exp(time)*np.sin(self.omega*time + x2/self.L)) + self.delta_alpha*x2*(self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L) + 20000*x1*x2*(x1 - 1)*(2*x1 - 1)*(x2 - 1)**2*np.cos(np.pi*time))*np.sin(self.omega*time + x2/self.L) + x2*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))*(40000*self.L**2*x1*(x1 - 1)*(2*x1 - 1)*(3*x2 - 2)*np.cos(np.pi*time) - self.delta_alpha*self.omega*np.exp(time)*np.cos(self.omega*time + x2/self.L))))*np.exp(-time)/(self.L**2*(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L))**3))

    def du21(self, time, x1, x2, x3):
        return (-20000*x2**2*(x2 - 1)**2*(2*x1*(x1 - 1) + x1*(2*x1 - 1) + (x1 - 1)*(2*x1 - 1))*np.exp(-time)*np.cos(np.pi*time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L)))

    def du211(self, time, x1, x2, x3):
        return (-120000*x2**2*(2*x1 - 1)*(x2 - 1)**2*np.exp(-time)*np.cos(np.pi*time)/(self.alpha0 + self.delta_alpha*np.sin(self.omega*time + x2/self.L)))