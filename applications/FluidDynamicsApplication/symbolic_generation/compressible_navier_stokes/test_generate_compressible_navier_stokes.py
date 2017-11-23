from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

from params_dict import params
import ConvectiveFlux
import DiffusiveFlux
import SourceTerm
import StabilizationMatrix

dim = params["dim"]  
dimes = dim+2 					        # Dimension of the vector of Unknowns
do_simplifications = False
#dim_to_compute = "Both"                # Spatial dimensions to compute. Options:  "2D","3D","Both"
mode = "c"                              # Output mode to a c++ file

if(dim == 2):
    nnodes = 3
elif(dim == 3):
    nnodes = 4

impose_partion_of_unity = False
N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)
  
# Unknown fields definition (Used later for the gauss point interpolation)
U = DefineMatrix('U',nnodes,dimes)	     # Vector of Unknowns ( Density,Velocity[dim],Total Energy )
Un = DefineMatrix('Un',nnodes,dimes)     # Vector of Unknowns one step back
Unn = DefineMatrix('Unn',nnodes,dimes)   # Vector of Unknowns two steps back
r = DefineVector('r',nnodes)             # Sink term

# Test functions defintiion
w = DefineMatrix('w',nnodes,dimes)	     # Variables field test

# External terms definition
f_ext = DefineMatrix('f_ext',nnodes,dim) # Forcing term 

# Definition of other symbols
bdf0 = Symbol('bdf0')                    # Backward differantiation coefficients
bdf1 = Symbol('bdf1')
bdf2 = Symbol('bdf2')


### Construction of the variational equation

Ug = DefineVector('Ug',dim+2)			# Dofs vector
H = DefineMatrix('H',dim+2,dim)			# Gradient of U
f = DefineVector('f',dim)			    # Body force vector
rg = Symbol('rg', positive = True)		# Source/Sink term
V = DefineVector('V',dim+2)			    # Test function
Q = DefineMatrix('Q',dim+2,dim)			# Gradient of V
acc = DefineVector('acc',dimes)         # Derivative of Dofs/Time

S = SourceTerm.computeS(f,rg,params)
#SourceTerm.printS(S,params)
A = ConvectiveFlux.computeA(Ug,params)
#ConvectiveFlux.printA(A,params)
K = DiffusiveFlux.computeK(Ug,params)
#DiffusiveFlux.printK(K,params)
Tau = StabilizationMatrix.computeTau(Ug,params)
#StabilizationMatrix.printTau(Tau,params)

print("\nCompute the residual\n")

## Nonlinear operator definition
L = DefineVector('L',dim+2)		       # Nonlinear operator
res = DefineVector('res',dim+2)		   # Residual definition

l1 = Matrix(zeros(dim+2,1))		       # Convective Matrix*Gradient of U
tmp = []
for j in range(0,dim):
    tmp = A[j]*H[:,j]
    l1 +=tmp

l2 = Matrix(zeros(dim+2,1))		       # Diffusive term

for s in range(0,dim+2):
    for k in range(0,dim):
        kinter = K[k]			       # Intermediate matrix
        for j in range(0,dim):
            ksmall = kinter[j]		   # Intermediate 2 matrix
            for l in range(0,dim+2):
                for m in range(0,dim+2):
                    for n in range(0,dim+2):
                        l2[s] += diff(ksmall[l,m],Ug[n])*H[n,k]*H[s,j]
                        #print("l",l,"m",m,"n",n,":\n",l2[s],"\n\n\n\n\n\n"

l3 = S*Ug				               # Source term
L = l1-l2-l3

## Redisual definition
res = -acc - L		
#res = -acc 

print("\nCompute the adjoint\n")
## Nonlinear adjoint operator definition

L_adj = DefineVector('L_adj',dim+2)	   # Nonlinear adjoint operator

m1 = Matrix(zeros(dim+2,1))		       # Convective term
for s in range(0,dim+2):
    for j in range(0,dim):
        A_T = A[j].transpose()
        for l in range(0,dim+2):
            for m in range(0,dim+2):
                m1[s] -= A_T[l,m]*Q[s,j]
                for n in range(0,dim+2):
                    m1[s] -= diff(A_T[l,m],Ug[n])*H[n,j]*V[s]
                  
m2 = Matrix(zeros(dim+2,1))		       # Diffusive term

for s in range(0,dim+2):
    for k in range(0,dim):
        kinter = K[k]
        for j in range(0,dim):
            ksmall = kinter[j].transpose()
            for l in range(0,dim+2):
                for m in range(0,dim+2):
                    for n in range(0,dim+2):
                        m2[s] -= diff(ksmall[l,m],Ug[n])*H[n,j]*Q[s,k]

m3 = -S.transpose()*V			        # Source term
L_adj = m1+m2+m3


L_adj[0] = 1;
L_adj[1] = 1;
L_adj[2] = 1;
L_adj[3] = 1;		  
## Variational Formulation - Final equation
print("\nCompute Variational Formulation - Final equation\n")

n1 = V.transpose()*acc		            # Mass term - FE scale

for i in range(0,dim):
    tmp += A[i]*H[:,i]
n2 = V.transpose()*tmp			       # Convective term - FE scale

tmp = Matrix(zeros(dim+2,1))
n3 = Matrix(zeros(1,1))			       # Diffusive term - FE scale
for k in range(0,dim):
    kinter= K[k]
    for j in range(0,dim):
        tmp += kinter[j]*H[:,j] 
    n3 += Q[:,k].transpose()*tmp

n4 = -V.transpose()*(S*Ug)		       # Source term - FE scale

n5 = L_adj.transpose()*(Tau*res)	   # VMS_adjoint - Subscales


rv = n1+n2+n3+n4+n5 			       # VARIATIONAL FORMULATION - FINAL EQUATION
#rv = n5 			       # VARIATIONAL FORMULATION - FINAL EQUATION

#print(rv)
### Substitution of the discretized values at the gauss points
print("\nSubstitution of the discretized values at the gauss points\n")

if(dim == 2 and nnodes==3): #I did it for 2D!
    U[0,0] = 1.177; U[0,1] = 1;    U[0,2] = -1; U[0,3] = 50;
    U[1,0] =  1.177; U[1,1] =  -1;    U[1,2] =  0; U[1,3] =  100;
    U[2,0] =  1.177; U[2,1] =  2;    U[2,2] =  1; U[2,3] =  0;
    Un = Matrix(zeros(nnodes,dimes));
    Unn = Matrix(zeros(nnodes,dimes));
    
    geom = [[0,0],[5,1],[3,5]]
    x10 = geom[1][0] - geom[0][0]
    y10 = geom[1][1] - geom[0][1]
    x20 = geom[2][0] - geom[0][0]
    y20 = geom[2][1] - geom[0][1]
    detJ = x10 * y20-y10 * x20;

    N[0] = 0.333333333333333;
    N[1] = 0.333333333333333;
    N[2] = 0.333333333333333;
    
    DN[0,0] = -y20 + y10; DN[0,1] = x20 - x10;
    DN[1,0] = y20; DN[1,1] = -x20;
    DN[2,0] =-y10; DN[2,1] = x10;
    
    for i in range(0,2):
        for j in range(0,1):
            DN[i,j]/=detJ
            
    f_ext[0,0] = 0; f_ext[0,1] = 0;
    f_ext[1,0] = 0; f_ext[1,1] = 1000;
    f_ext[2,0] = 0; f_ext[2,1] = 0;
    
    # Definition of other symbols
    bdf0 = 2/3                    # Backward differantiation coefficients
    bdf1 = 4/3
    bdf2 = -1/3
    r[0] = 33; r[1] = -33; r[2] = 33; 
    
    #w[0,0] = 300; w[0,1] = 1;    w[0,2] = -1; w[0,3] = 50;
    #w[1,0] =  300; w[1,1] =  -1;    w[1,2] =  0; w[1,3] =  100;
    #w[2,0] =  300; w[2,1] =  2;    w[2,2] =  1; w[2,3] =  0;

## Data interpolation at the gauss points
U_gauss = U.transpose()*N
#print('\nU_gauss = ', U_gauss)

w_gauss = w.transpose()*N
f_gauss = f_ext.transpose()*N
acc_gauss = (bdf0*U+bdf1*Un+bdf2*Unn).transpose()*N
r_gauss = (r.transpose()*N)[0]      

## Gradients computation
grad_U = DfjDxi(DN,U).transpose()
grad_w = DfjDxi(DN,w).transpose()

SubstituteMatrixValue(rv, Ug, U_gauss)
SubstituteMatrixValue(rv, acc, acc_gauss)
SubstituteMatrixValue(rv, H, grad_U)
SubstituteMatrixValue(rv, V, w_gauss)
SubstituteMatrixValue(rv, Q, grad_w)
SubstituteMatrixValue(rv, f, f_gauss)
SubstituteScalarValue(rv, rg, r_gauss)

print("\n\n\n\n")
print(rv)

print("\n\n\n\n")


## Compute LHS and RHS

rhs = Compute_RHS(rv.copy(), w, do_simplifications)
rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)
print("\nRHS:\n")
print("\n0:",rhs[0])
print("\n1:",rhs[1])
print("\n2:",rhs[2])

lhs = Compute_LHS(rhs, w, U, do_simplifications) # Compute the LHS
lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

print("\nLHS:\n")
for i in range(w.shape[0]):
    for j in range (U.shape[0]):
        print("\n",i,j,"=",lhs[i,j])
        


