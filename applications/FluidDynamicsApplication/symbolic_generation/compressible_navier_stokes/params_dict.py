## Dictionary of the constant parameters used in the Variational Formulation
from sympy import *

params = {
            "dim": 3,			                            # Dimension
            "mu": Symbol('mu', positive = True),			# Dynamic viscosity 
            "nu" :Symbol('nu', positive = True),			# Kinematic viscosity (mu/rho)
            "h" : Symbol('h', positive = True),	        	# Element size
            "lambda" : Symbol('l', positive = True),         # Thermal Conductivity of the fluid
            "c_v" : Symbol('cv', positive = True),			# Specific Heat at Constant volume
            "c_p" : Symbol('cp', positive = True),			# Specific Heat at Constant Pressure
            "gamma": Symbol('y',positive = True),			# Gamma (Cp/Cv) 
            "c_1" : Symbol('stab_c1', positive = True),			# Algorithm constant
            "c_2" : Symbol('stab_c2', positive = True),			# Algorithm constant
            
	}
#print (params)
'''
## 2D test
params = {
            "dim": 2,			                            # Dimension
            "mu": 0.798,			# Dynamic viscosity 
            "nu" :0.801,			# Kinematic viscosity (mu/rho)
            "e" : 125.75,	        	# Enthalpy
            "lambda" : 0.591,         # Thermal Conductivity of the fluid
            "c_v" : 74.53,			# Specific Heat at Constant volume
            "c_p" : 4.1796,			# Specific Heat at Constant Pressure
            "gamma": 0.056079,			# Gamma (Cp/Cv) 
            "c_1" : 4,			# Algorithm constant
            "c_2" : 2,			# Algorithm constant
            
	}
#print (params)
'''
