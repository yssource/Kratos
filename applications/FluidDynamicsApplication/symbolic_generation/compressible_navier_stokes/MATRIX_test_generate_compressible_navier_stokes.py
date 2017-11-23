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
dim_to_compute = "2D"                # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
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
    Nval,DNval = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)
    
    # Unknown fields definition (Used later for the gauss point interpolation)
    U = DefineMatrix('U',nnodes,dimes, real=True)	     # Vector of Unknowns ( Density,Velocity[dim],Total Energy )
    Uval = DefineMatrix('U',nnodes,dimes, real=True)	 
    Un = DefineMatrix('Un',nnodes,dimes, real=True)     # Vector of Unknowns one step back
    Unn = DefineMatrix('Unn',nnodes,dimes, real=True)   # Vector of Unknowns two steps back
    r = DefineVector('r',nnodes, real=True)             # Sink term

    # Test functions defintiion
    w = DefineMatrix('w',nnodes,dimes, real=True)	     # Variables field test

    # External terms definition
    f_ext = DefineMatrix('f_ext',nnodes,dim, real=True) # Forcing term 
    f_ext_val = DefineMatrix('f_ext',nnodes,dim, real=True) # Forcing term 

    # Definition of other symbols
    bdf0 = Symbol('bdf0', real=True)                    # Backward differantiation coefficients
    bdf1 = Symbol('bdf1', real=True)
    bdf2 = Symbol('bdf2', real=True)

    ### Construction of the variational equation
    Ug = DefineVector('Ug',dim+2, real=True)			# Dofs vector
    H = DefineMatrix('H',dim+2,dim, real=True)			# Gradient of U
    f = DefineVector('f',dim, real=True)			    # Body force vector
    rg = Symbol('rg', real=True)		# Source/Sink term
    V = DefineVector('V',dim+2, real=True)			    # Test function
    Q = DefineMatrix('Q',dim+2,dim, real=True)			# Gradient of V
    acc = DefineVector('acc',dimes, real=True)         # Derivative of Dofs/Time

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
    L = DefineVector('L',dim+2, real=True)		       # Nonlinear operator
    res = DefineVector('res',dim+2, real=True)		   # Residual definition

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

    L_adj = DefineVector('L_adj',dim+2, real=True)	   # Nonlinear adjoint operator

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

    ## Data interpolation at the gauss points
    U_gauss = U.transpose()*N
    w_gauss = w.transpose()*N
    f_gauss = f_ext.transpose()*N
    acc_gauss = (bdf0*U+bdf1*Un+bdf2*Unn).transpose()*N
    r_gauss = (r.transpose()*N)[0]    
    #print(r_gauss)

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
    #print(rv)
   
    dofs = Matrix(zeros(nnodes*(dim+2),1), real=True)
    testfunc = Matrix(zeros(nnodes*(dim+2),1), real=True)
    
    for i in range(0,nnodes):
         for j in range(0,dim+2):
            dofs[i*(dim+2)+j] = U[i,j]
            testfunc[i*(dim+2)+j] = w[i,j]
    
    ## Compute LHS and RHS
    print("\nCompute RHS\n")
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    #print(rhs)
    
    #rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)
    
    '''
    ## TEST TO CHECK NUMERICAL RESULT IN RHS
    print("\nNUMERICAL RHS:\n")
    dofsval = Matrix(zeros(nnodes*(dim+2),1))
    testfuncval = Matrix(zeros(nnodes*(dim+2),1))
    
    if(dim == 2 and nnodes==3): #I did it for 2D!
        dofsval[0] = 1.177; dofsval[1] = 1;    dofsval[2] = -1; dofsval[3] = 50;
        dofsval[4] =  1.177; dofsval[5] =  -1;    dofsval[6] =  0; dofsval[7] =  100;
        dofsval[8] =  1.177; dofsval[9] =  2;    dofsval[10] =  1; dofsval[11] =  0;
        Unval = Matrix(zeros(nnodes,dimes));
        Unnval = Matrix(zeros(nnodes,dimes));
        rval = Matrix(zeros(nnodes,1));
        
        testfuncval[0] = 300; testfuncval[1] = 1;    
        testfuncval[2] = -1; testfuncval[3] = 50;
        testfuncval[4] =  300; testfuncval[5] =  -1;   testfuncval[6] =  0; testfuncval[7] =  100;
        testfuncval[8] =  300; testfuncval[8] =  2;    testfuncval[10] =  1; testfuncval[11] =  0;
        
        geom = [[0,0],[5,1],[3,5]]
        x10 = geom[1][0] - geom[0][0]
        y10 = geom[1][1] - geom[0][1]
        x20 = geom[2][0] - geom[0][0]
        y20 = geom[2][1] - geom[0][1]
        detJ = x10 * y20-y10 * x20;

        Nval[0] = 0.33;
        Nval[1] = 0.33;
        Nval[2] = 0.33;
        
        DNval[0,0] = -y20 + y10; DNval[0,1] = x20 - x10;
        DNval[1,0] = y20; DNval[1,1] = -x20;
        DNval[2,0] =-y10; DNval[2,1] = x10;
        
        for i in range(0,2):
            for j in range(0,1):
                DNval[i,j]/=detJ
                
        f_ext_val[0,0] = 0; f_ext_val[0,1] = 0;
        f_ext_val[1,0] = 0; f_ext_val[1,1] = 1000;
        f_ext_val[2,0] = 0; f_ext_val[2,1] = 0;
        
        # Definition of other symbols
        bdf0val = 2/3                    # Backward differantiation coefficients
        bdf1val = 4/3
        bdf2val = -1/3
        rval[0] = 33; rval[1] = 22; rval[2] = -12; 
        
        
        print("\nSubstitute testfunc:\n")
        SubstituteMatrixValue(rhs, testfunc, testfuncval)
        for i in range(0,rhs.shape[0]):
            print("rhs(",i,"):=",rhs[i],"\t\t")
            print("\n")
        print("\nSubstitute Un, Unn:\n")
        SubstituteMatrixValue(rhs, Un, Unval)
        SubstituteMatrixValue(rhs, Unn, Unnval)
        for i in range(0,rhs.shape[0]):
            print("rhs(",i,"):=",rhs[i],"\t\t")
            print("\n")
        print("\nSubstitute N,Dn:\n")
        SubstituteMatrixValue(rhs, N, Nval)
        SubstituteMatrixValue(rhs, DN, DNval)
        for i in range(0,rhs.shape[0]):
            print("rhs(",i,"):=",rhs[i],"\t\t")
            print("\n")
        print("\nSubstitute bdf,rval:\n")
        SubstituteScalarValue(rhs, bdf0, bdf0val)
        SubstituteScalarValue(rhs, bdf1, bdf1val)
        SubstituteScalarValue(rhs, bdf2, bdf2val)
        SubstituteMatrixValue(rhs, r, rval)
        for i in range(0,rhs.shape[0]):
            print("rhs(",i,"):=",rhs[i],"\t\t")
            print("\n")
        print("\nSubstitute fext:\n")
        SubstituteMatrixValue(rhs, f_ext, f_ext_val)
        print("\nSubstitute DOFS:\n")
        SubstituteMatrixValue(rhs, dofs, dofsval)
        
        print("\nNUMERICAL RHS:\n")
        for i in range(0,rhs.shape[0]):
            print("rhs(",i,"):=",rhs[i],"\t\t")
            print("\n")
    '''
        
    print("\nCompute LHS\n")
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)
   
    ## TEST TO CHECK NUMERICAL RESULT IN LHS
    print("\nLHS:\n")
    dofsval = Matrix(zeros(nnodes*(dim+2),1), real=True)
    testfuncval = Matrix(zeros(nnodes*(dim+2),1), real=True)
    
    if(dim == 2 and nnodes==3): #I did it for 2D!
        dofsval[0] = 1.177; dofsval[1] = 1;    dofsval[2] = -1; dofsval[3] = 50;
        dofsval[4] =  1.177; dofsval[5] =  -1;    dofsval[6] =  0; dofsval[7] =  100;
        dofsval[8] =  1.177; dofsval[9] =  2;    dofsval[10] =  1; dofsval[11] =  0;
        Unval = Matrix(zeros(nnodes,dimes));
        Unnval = Matrix(zeros(nnodes,dimes));
        rval = Matrix(zeros(nnodes,1));
        
        testfuncval[0] = 300; testfuncval[1] = 1;    
        testfuncval[2] = -1; testfuncval[3] = 50;
        testfuncval[4] =  300; testfuncval[5] =  -1;   testfuncval[6] =  0; testfuncval[7] =  100;
        testfuncval[8] =  300; testfuncval[8] =  2;    testfuncval[10] =  1; testfuncval[11] =  0;
        
        geom = [[0,0],[5,1],[3,5]]
        x10 = geom[1][0] - geom[0][0]
        y10 = geom[1][1] - geom[0][1]
        x20 = geom[2][0] - geom[0][0]
        y20 = geom[2][1] - geom[0][1]
        detJ = x10 * y20-y10 * x20;

        Nval[0] = 0.33;
        Nval[1] = 0.33;
        Nval[2] = 0.33;
        
        DNval[0,0] = -y20 + y10; DNval[0,1] = x20 - x10;
        DNval[1,0] = y20; DNval[1,1] = -x20;
        DNval[2,0] =-y10; DNval[2,1] = x10;
        
        for i in range(0,2):
            for j in range(0,1):
                DNval[i,j]/=detJ
                
        f_ext_val[0,0] = 0; f_ext_val[0,1] = 0;
        f_ext_val[1,0] = 0; f_ext_val[1,1] = 1000;
        f_ext_val[2,0] = 0; f_ext_val[2,1] = 0;
        
        # Definition of other symbols
        bdf0val = 2/3                    # Backward differantiation coefficients
        bdf1val = 4/3
        bdf2val = -1/3
        rval[0] = 33; rval[1] = 22; rval[2] = -12; 
            
        print("\nSubstitute DOFS:\n")
        SubstituteMatrixValue(lhs, dofs, dofsval)
        print("\nSubstitute N,Dn:\n")
        SubstituteMatrixValue(lhs, N, Nval)
        SubstituteMatrixValue(lhs, DN, DNval) 
        for i in range(0,nnodes*(dim+2)):
            for j in range(0,nnodes*(dim+2)):
                print("lhs(",i,",",j,"):=",lhs[i,j],"\t\t")
            print("\n")
        print("\nSubstitute testfunc:\n")
        SubstituteMatrixValue(lhs, testfunc, testfuncval)
        print("\nSubstitute Un, Unn:\n")
        SubstituteMatrixValue(lhs, Un, Unval)
        SubstituteMatrixValue(lhs, Unn, Unnval)
        print("\nSubstitute bdf,rval:\n")
        SubstituteScalarValue(lhs, bdf0, bdf0val)
        SubstituteScalarValue(lhs, bdf1, bdf1val)
        SubstituteScalarValue(lhs, bdf2, bdf2val)
        SubstituteMatrixValue(lhs, r, rval)
        print("\nSubstitute fext:\n")
        SubstituteMatrixValue(lhs, f_ext, f_ext_val)
        
        
        print("\nNUMERICAL LHS:\n")
        for i in range(0,nnodes*(dim+2)):
            for j in range(0,nnodes*(dim+2)):
                print("lhs(",i,",",j,"):=",lhs[i,j],"\t\t")
            print("\n") 

'''
    ## READING TEMPLATE FILE
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
'''



