from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

def computeTau(dofs,params):
    print("\nCompute Stabilization Matrix\n")
    dim = params["dim"]				# spatial dimensions
    ## Unknown field definition
    Tau = DefineMatrix('Tau',dim+2,dim+2)	# Stabilization matrix 
    
    ## Data interpolation to the Gauss points
    Ug = dofs

    ## Other symbols definitions
    y = params["gamma"]				# Gamma (Cp/Cv)
    Cp = params["c_p"]				# Specific Heat at Constant Pressure
    nu = params["nu"]				# Kinematic viscosity (mu/rho)
    l = params["lambda"]			# Thermal Conductivity of the fluid
    h = params["h"]				# Element size
    c1 = params["stab_c1"]				# Algorithm constant
    c2 = params["stab_c2"]				# Algorithm constant
    
    ## c - Speed of Sound definition
    c_tmp = Ug[dim+1]/Ug[0]
    for i in range (0,dim):
        c_tmp += -Ug[i+1]**2/(2*Ug[0]**2)
    c_g = sqrt((y*(y-1)*c_tmp))
 
    ## Tau - Stabilization Matrix definition
    Tau = zeros(dim+2,dim+2)
    
    tau1 = 0.0
    for i in range(0,dim):
        tau1 +=(Ug[i+1]/Ug[0])**2
    tau1 = (sqrt(tau1)+c_g)*c2/h
    tau2 = c1*nu/h**2+tau1
    tau3 = c1*l/(Ug[0]*Cp*h**2)+tau1
    
    Tau[0,0] = tau1
    for i in range (0,dim):
    	Tau[i+1,i+1] = tau2
    Tau[dim+1,dim+1] = tau3
    
    Tau = Tau.inv()
    '''
    for i in range (0, dim+2):
        for j in range (0,dim+2):
            print(i,j,"=",Tau[i,j],"\n")
    '''
    return(Tau)

def computeSymbolicTau(dofs,params):
    print("\nCompute Stabilization Matrix\n")
    dim = params["dim"]				# spatial dimensions
    ## Unknown field definition
    if(dim == 2):
        nnodes = 3
    elif(dim == 3):
        nnodes = 4

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)
   
    Tau = zeros(dim+2,dim+2)
    tau1 = Symbol('tau1')
    tau2 = Symbol('tau2')
    tau3 = Symbol('tau3')
    

    Tau[0,0] = tau1
    for i in range (0,dim):
        Tau[i+1,i+1] = tau2
    Tau[dim+1,dim+1] = tau3

    Tau = Tau.inv()
    return(Tau)
    
def printTau(Tau, params):
    dim = params["dim"]				#spatial dimensions
    print("The Stabilization term matrix is:\n")
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("Tau[",i,",",j,"]=",Tau[i,j],"\n")

    return 0
    
    
    
    
    
    
    
    
    
    
