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
BlockSize = dim+2 					        # Dimension of the vector of Unknowns
do_simplifications = False
dim_to_compute = "2D"                # Spatial dimensions to compute. Options:  "2D","3D","Both"
mode = "c"                              # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open("compressible_navier_stokes_cpp_template.cpp")
outstring = templatefile.read()

for dim in dim_vector:

    if(dim == 2):
        nnodes = 3
    elif(dim == 3):
        nnodes = 4
    
    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)
    
    # Unknown fields definition (Used later for the gauss point interpolation)
    U = DefineMatrix('U',nnodes,BlockSize)	     # Vector of Unknowns ( Density,Velocity[dim],Total Energy )
    Un = DefineMatrix('Un',nnodes,BlockSize)     # Vector of Unknowns one step back
    Unn = DefineMatrix('Unn',nnodes,BlockSize)   # Vector of Unknowns two steps back
    r = DefineVector('r',nnodes)             # Sink term    #COMMENT for manufactured solution

    # Test functions defintiion
    w = DefineMatrix('w',nnodes,BlockSize)	     # Variables field test

    # External terms definition
    f_ext = DefineMatrix('f_ext',nnodes,dim) # Forcing term #COMMENT for manufactured solution

    # Definition of other symbols
    bdf0 = Symbol('bdf0')                    # Backward differantiation coefficients
    bdf1 = Symbol('bdf1')
    bdf2 = Symbol('bdf2')


    ### Construction of the variational equation

    Ug = DefineVector('Ug',BlockSize)			# Dofs vector
    H = DefineMatrix('H',BlockSize,dim)			# Gradient of U
    f = DefineVector('f',dim)			    # Body force vector
    rg = Symbol('rg', positive = True)		# Source/Sink term
    V = DefineVector('V',BlockSize)			    # Test function
    Q = DefineMatrix('Q',BlockSize,dim)			# Gradient of V
    acc = DefineVector('acc',BlockSize)         # Derivative of Dofs/Time
    G = DefineMatrix('G',BlockSize,dim)		# Diffusive Flux matrix
    
    S = SourceTerm.computeS(f,rg,params)
    #SourceTerm.printS(S,params)
    A = ConvectiveFlux.computeA(Ug,params)
    #ConvectiveFlux.printA(A,params)
    G = DiffusiveFlux.computeK(Ug,params,H,G)
    #DiffusiveFlux.printK(G,params)
    Tau = StabilizationMatrix.computeTau(params)
    #StabilizationMatrix.printTau(Tau,params)
    
    ## Nonlinear operator definition   
    l1 = Matrix(zeros(dim+2,1))		       # Convective Matrix*Gradient of U
    A_small = []
    for j in range(0,dim):
        A_small = A[j]
        for ll in range(BlockSize):
            for mm in range(BlockSize):
                l1[ll] += A_small[ll,mm]*H[mm,j]
    '''
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
    
    '''
    l3 = S*Ug				               # Source term
    print("\nCompute Non-linear operator\n")
    L = l1-l3 #-l2                          # Nonlinear operator

    ## Residual definition     
    res = -acc - L		
   
    ## Nonlinear adjoint operator definition  
    m1 = Matrix(zeros(dim+2,1))		       # Convective term
    psi = Matrix(zeros(dim+2,dim))

    for j in range(0,dim):
        A_T = A[j].transpose()
        for l in range(0,dim+2):
            for m in range(0,dim+2):
                psi[l,j] += A_T[l,m]*Q[m,j]                 
                for n in range(0,dim+2):
                    psi[l,j] +=diff(A_T[l,m],Ug[n])*H[n,j]*V[m]   

        #print(psi)

    for s in range(0,dim+2):
        for j in range(0,dim):
            m1[s] += psi[s,j]

    '''
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
    '''
    m3 = S.transpose()*V			        # Source term
    L_adj = -m1-m3 #+m2                      # Nonlinear adjoint operator
         
    ## Variational Formulation - Final equation
    n1 = V.transpose()*acc		            # Mass term - FE scale
     
    temp = zeros(dim+2,1)
    A_smalll = []
    for i in range(0,dim):
        A_smalll = A[i]
        for ll in range(BlockSize):
            for mm in range(BlockSize):
                temp[ll] += A_smalll[ll,mm]*H[mm,i]
    
    n2 = V.transpose()*temp			       # Convective term - FE scale

    '''
    tmp = Matrix(zeros(dim+2,1))
    n3 = Matrix(zeros(1,1))			       
    
    for k in range(0,dim):
        kinter= K[k]
        for j in range(0,dim):
            tmp += kinter[j]*H[:,j] 
        n3 += Q[:,k].transpose()*tmp
    '''

    n3 = Matrix(zeros(1,1))                 # Diffusive term - FE scale
    for j in range(0,dim):
        for k in range(BlockSize):
            n3[0,0] += Q[k,j]*(-G[k,j])

    n4 = -V.transpose()*(S*Ug)		       # Source term - FE scale
    
    n5 = L_adj.transpose()*(Tau*res)	   # VMS_adjoint - Subscales 
 
    print("\nCompute Variational Formulation\n")
    rv = n1+n2+n3+n4+n5 			       # VARIATIONAL FORMULATION - FINAL EQUATION
  

    ### Substitution of the discretized values at the gauss points
    print("\nSubstitution of the discretized values at the gauss points\n")
    
    ## Data interpolation at the gauss points
    U_gauss = U.transpose()*N
    w_gauss = w.transpose()*N
    f_gauss = f_ext.transpose()*N                     #COMMENT for manufactured solution
    acc_gauss = (bdf0*U+bdf1*Un+bdf2*Unn).transpose()*N
    r_gauss = (r.transpose()*N)[0]                   #COMMENT for manufactured solution   
    #r_gauss = Symbol('r_gauss', positive = True)     #USED FOR MANUFACTURED SOLUTION

    ## Gradients computation
    grad_U = DfjDxi(DN,U).transpose()
    grad_w = DfjDxi(DN,w).transpose()
    
    SubstituteMatrixValue(rv, Ug, U_gauss)
    SubstituteMatrixValue(rv, acc, acc_gauss)
    SubstituteMatrixValue(rv, H, grad_U)
    SubstituteMatrixValue(rv, V, w_gauss)
    SubstituteMatrixValue(rv, Q, grad_w)
    SubstituteMatrixValue(rv, f, f_gauss)       #COMMENT for manufactured solution
    SubstituteScalarValue(rv, rg, r_gauss)      #COMMENT for manufactured solution
    
    dofs = Matrix(zeros(nnodes*(dim+2),1))
    testfunc = Matrix(zeros(nnodes*(dim+2),1))
    for i in range(0,nnodes):
         for j in range(0,dim+2):
            dofs[i*(dim+2)+j] = U[i,j]
            testfunc[i*(dim+2)+j] = w[i,j]

    ## Compute LHS and RHS
    print("\nCompute RHS\n")
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)
        
    print("\nCompute LHS\n")
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)
    
    ## Reading Template File
    print("\nReading compressible_navier_stokes_cpp_template.cpp\n")
    templatefile = open("compressible_navier_stokes_cpp_template.cpp")
    outstring=templatefile.read()

    if(dim == 2):
            outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
            outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
            outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
            outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

## Write the modified template
print("\nWriting compressible_navier_stokes.cpp\n")
out = open("compressible_navier_stokes.cpp",'w')
out.write(outstring)
out.close()
