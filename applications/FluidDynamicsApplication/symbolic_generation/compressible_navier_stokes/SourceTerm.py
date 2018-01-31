from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

## Computation of the Source Matrix
def computeS(force,source,params):
    print("\nCompute Source Matrix \n")
    dim = params["dim"]				# spatial dimensions
    
    ## Unknown field definition
    #S = DefineMatrix('S',dim+2,dim+2)		# Reactive matrix (Source terms)
    f = force					            # Body force vector
    r = source			 		            # Heat Source/Sink Term  #Symbol('r', positive = True) 

    ## S - Reactive Matrix definition
    S = zeros(dim+2,dim+2)

    #0  0  0 0
    #fx 0  0 0
    #fy 0  0 0
    #r fx fy 0
    for i in range(1,dim+1):
        S[i,0] = f[i-1]
    S[dim+1,0] = r

    for j in range(1,dim+1):
        S[dim+1,j] = f[j-1]
        
    return S

## Printing of the Source Matrix   
def printS(S,params):
    dim = params["dim"]				#spatial dimensions
    print("The source term matrix is:\n")
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("S[",i,",",j,"]=",S[i,j],"\n")

    return 0
