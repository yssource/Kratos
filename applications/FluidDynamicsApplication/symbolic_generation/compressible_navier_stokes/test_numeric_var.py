from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

## TWO DIMENSIONAL TEST
# using one gauss point

#N = []
#DN_DX = []
#DN_DX = Matrix(zeros(3,2))	
#geom = DefineMatrix('geom',3,2)
geom = [[0,0],[5,1],[3,5]]

x10 = geom[1][0] - geom[0][0]
y10 = geom[1][1] - geom[0][1]
x20 = geom[2][0] - geom[0][0]
y20 = geom[2][1] - geom[0][1]

detJ = x10 * y20-y10 * x20;

DN_DX = [[-y20 + y10,x20 - x10],[y20, -x20],[-y10,x10]]

for i in range(0,2):
    for j in range(0,1):
        DN_DX[i][j]/=detJ

N = [[0.333333333333333],[0.333333333333333],[0.333333333333333]]
N[0] = 0.333333333333333;
N[1] = 0.333333333333333;
N[2] = 0.333333333333333;

Area = 0.5*detJ

def set_N():
    return(N)

def set_DN():
    return(DN_DX)
