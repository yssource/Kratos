from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

from params_dict import params

#def computeTau(dofs,params):
print("\nCompute Stabilization Matrix\n")
dim = params["dim"]				# spatial dimensions
## Unknown field definition
if(dim == 2):
    nnodes = 3
elif(dim == 3):
    nnodes = 4

impose_partion_of_unity = False
N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)
Tau = DefineMatrix('Tau',dim+2,dim+2)	# Stabilization matrix 

## Data interpolation to the Gauss points
Ug = DefineVector('Ug',dim+2)			# Dofs vector

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
c_g = (y*(y-1)*c_tmp)
print("\nSound sqrt\n\n",c_g)
#c_g = real_root((y*(y-1)*c_tmp))
#print("\nSound real_root\n",c_g)



## Tau - Stabilization Matrix definition
tmp = 0

for i in range (0,dim):  
    tmp += Ug[i+1]
    
Tau = zeros(dim+2,dim+2)

tau1 = c2*(abs(tmp/Ug[0])+c_g)/h
tau2 = c1*nu/h**2+tau1
tau3 = c1*l/(Ug[0]*Cp*h**2)+tau1

Tau[0,0] = tau1
for i in range (0,dim):
    Tau[i+1,i+1] = tau2
Tau[dim+1,dim+1] = tau3

Tau = Tau.inv()

U = DefineMatrix('U',nnodes,dim+2)
if(dim == 2 and nnodes==3): #I did it for 2D!
    U[0,0] = 1.2041	; U[0,1] = 1;    U[0,2] = -1; U[0,3] = 50;
    U[1,0] =  1.2041	; U[1,1] =  -1;    U[1,2] =  0; U[1,3] =  100;
    U[2,0] =  1.2041	; U[2,1] =  2;    U[2,2] =  1; U[2,3] =  0;
    
    geom = [[0,0],[5,1],[3,5]]
    x10 = geom[1][0] - geom[0][0]
    y10 = geom[1][1] - geom[0][1]
    x20 = geom[2][0] - geom[0][0]
    y20 = geom[2][1] - geom[0][1]
    detJ = x10 * y20-y10 * x20;

    N[0] = 0.333333333333333;
    N[1] = 0.333333333333333;
    N[2] = 0.333333333333333;

U_gauss =  U.transpose()*N
SubstituteScalarValue(c_g,"Ug_0",'mamma')
print(c_g)
c_g = sqrt(c_g)
print(c_g)
'''
SubstituteScalarValue(c_g,Ug,Ug_1,U_gauss[1])
SubstituteScalarValue(c_g,Ug,Ug_0,U_gauss[0])
SubstituteScalarValue(c_g,Ug,Ug_0,U_gauss[0])
SubstituteScalarValue(Tau, Ug, Ug_0,U_gauss[0])
'''
print(c_g)
'''
for i in range (0, dim+2):
    for j in range (0,dim+2):
        print(i,j,"=",Tau[i,j],"\n")
'''
    #return(Tau)
    
def printTau(Tau, params):
    dim = params["dim"]				#spatial dimensions
    print("The Stabilization term matrix is:\n")
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("Tau[",i,",",j,"]=",Tau[i,j],"\n")

    return 0
    
    
    
    
    
    
    
    
    
    
