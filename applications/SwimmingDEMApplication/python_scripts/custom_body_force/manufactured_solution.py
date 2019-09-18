import KratosMultiphysics

def CreateManufacturedSolution(custom_settings):
    return ManufacturedSolution(custom_settings)

class ManufacturedSolution(object):
    def __init__(self, settings):
        '''
        This is a base class to build manufactured fluid solutions.
        At least, it should return the body force and the velocity.
        The input viscosity is the DYNAMIC viscosity
        NOTE: the operators are implemented for the 2D case. It could be extended to the 3D case.
        '''

        default_settings = KratosMultiphysics.Parameters("""
            {
                "viscosity"   : 1.0e-2,
                "density"     : 1.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.rho = settings["density"].GetDouble()
        self.nu  = settings["viscosity"].GetDouble() / self.rho

    # Public methods

    def BodyForce(self, t, *x):
        return [self.body_force1(t, *x), self.body_force2(t, *x), self.body_force3(t, *x)]

    def Velocity(self, t, *x):
        return [self.u1(t, *x), self.u2(t, *x), self.u3(t, *x)]

    def Pressure(self, t, *x):
        return self.p(t, *x)

    # Operators
    def body_force1(self, t, *x):
        return self.du1dt(t, *x) + self.convective1(t, *x) + 1 / self.rho * self.press_grad1(t, *x) - self.nu * self.laplacian1(t, *x)

    def body_force2(self, t, *x):
        return self.du2dt(t, *x) + self.convective2(t, *x) + 1 / self.rho * self.press_grad2(t, *x) - self.nu * self.laplacian2(t, *x)

    def body_force3(self, t, *x):
        return 0.0

    def convective1(self, t, *x):
        return  self.u1(t, *x) * self.du11(t, *x) + self.u2(t, *x) * self.du12(t, *x)

    def convective2(self, t, *x):
        return  self.u1(t, *x) * self.du21(t, *x) + self.u2(t, *x) * self.du22(t, *x)

    def laplacian1(self, t, *x):
        return  self.du111(t, *x) + self.du122(t, *x)

    def laplacian2(self, t, *x):
        return  self.du211(t, *x) + self.du222(t, *x)

    def press_grad1(self, t, *x):
        return self.dp1(t, *x)

    def press_grad2(self, t, *x):
        return self.dp2(t, *x)

    # Velocity and derivatives

    def u1(self, t, *x):
        """ Velocity
        """
        raise Exception("Method not implemented")

    def u2(self, t, *x):
        """ Velocity
        """
        raise Exception("Method not implemented")

    def u3(self, t, *x):
        return 0.0

    def du1dt(self, t, *x):
        raise Exception("Method not implemented")

    def du2dt(self, t, *x):
        raise Exception("Method not implemented")

    def du11(self, t, *x):
        raise Exception("Method not implemented")

    def du12(self, t, *x):
        raise Exception("Method not implemented")

    def du21(self, t, *x):
        raise Exception("Method not implemented")

    def du22(self, t, *x):
        raise Exception("Method not implemented")

    def du111(self, t, *x):
        raise Exception("Method not implemented")

    def du122(self, t, *x):
        raise Exception("Method not implemented")

    def du211(self, t, *x):
        raise Exception("Method not implemented")

    def du222(self, t, *x):
        raise Exception("Method not implemented")

    # Pressure and derivatives

    def p(self, t, *x):
        '''
        By default, pressure is 0
        '''
        return 0.0

    def dp1(self, t, *x):
        '''
        By default, pressure is 0
        '''
        return 0.0

    def dp2(self, t, *x):
        '''
        By default, pressure is 0
        '''
        return 0.0

    #Fluid fraction

    def alpha(self, t, *x):
        '''
        By default, alpha is 1
        '''
        return 1.0