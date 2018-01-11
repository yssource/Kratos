from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

## Computation of the Diffusive Matrix
def computeK(dofs,params):
    print("\nCompute Diffusive Matrix\n")
    dim = params["dim"]				# spatial dimensions
    ## Unknown fields definition
    H = DefineMatrix('H',dim+2,dim)		# Gradient of U
    G = DefineMatrix('G',dim+2,dim)		# Diffusive Flux matrix 
    tau = DefineMatrix('tau',dim,dim)		# Shear stress tensor for Newtonian fluid
    q = DefineVector('q',dim)			# Heat flux vector
    
    ## Other simbols definition
    Cv = params["c_v"]				# Specific Heat at Constant volume
    y = params["gamma"]				# Gamma (Cp/Cv) 
    mu  = params["mu"]         			# Dynamic viscosity 
    l = params["lambda"]			# Thermal Conductivity of the fluid
        
    ## Data interpolation to the Gauss points
    Ug = dofs

    ## Pgauss - Pressure definition
    pg = (y-1)*Ug[dim+1]
    for i in range(0,dim):
        pg += (y-1)*(-Ug[i+1]*Ug[i+1]/(2*Ug[0]))
    
    ## Tau - Shear stress tensor definition

    for i in range(0,dim):
        for j in range(i,dim):#NOT SURE
            if i!=j:
               tau[i,j] = mu/Ug[0]*(H[i+1,j]+H[j+1,i])-mu/Ug[0]**2*(Ug[i+1]*H[0,j]+Ug[j+1]*H[0,i])
            if i==j:
               tau[i,j]= 2*mu/Ug[0]*H[i+1,i]-2*mu/Ug[0]**2*Ug[i+1]*H[0,i]
               for k in range(0,dim):
                   tau[i,j]+= -2*mu/(3*Ug[0])*H[k+1,k]+2*mu/(3*Ug[0]**2)*Ug[k+1]*H[0,k]
    
    for i in range(1,dim):
        for j in range(0,dim-1):
            if j!=i:
               tau[i,j] = tau[j,i]
               
    ## q - Heat flux vector definition
    for i in range(0,dim):
        q[i] = l*Ug[dim+1]/(Ug[0]**2*Cv)*H[0,i]-(l*H[dim+1,i])/(Ug[0]*Cv)
        for j in range(0,dim):
            q[i] += -l*Ug[j+1]**2/(Cv*Ug[0]**3)*H[0,i]+l/(Ug[0]**2*Cv)*Ug[j+1]*H[j+1,i] 
    #NB!!!There is an error in the definition of q[i] in the research proposal. The second term of the equation has an opposite sign!!!NB#
       
    ## G - Diffusive Matrix definition 
    for j in range(0,dim):
        G[0,j]= 0 			#Mass equation related
       
    for i in range(1,dim+1):
        for j in range(0,dim):
            G[i,j]=-tau[i-1,j]		#Moment equation related
    
    for j in range(0,dim):
        G[dim+1,j] = q[j]
        for k in range(0,dim):
            G[dim+1,j] += -Ug[k+1]*tau[k,j]/Ug[0]
    
    ## K - Jacobian Diffusive Matrix definition
    K = []				#Final 5*5*3 tensor 			
    # k:index of H(moving over colomns), j:index of G(moving over colomns)	
      
    for k in range(0,dim):
        ksmall = []			#Intermediate 5*5*3 tensor
        for j in range(0,dim):
            tmp = DefineMatrix('tmp',dim+2,dim+2)
            for l in range(0,dim+2):
                for m in range(0,dim+2):
                    #tmp[l,m] = diff(-G[l,j],H[m,k])
                    tmp[l,m] = diff(G[l,k],H[m,j])
                    #print("\n\n",l,m,"=",tmp[l,m],"\n\n")
            
            ksmall.append(tmp)
        
        K.append(ksmall)     
    return K
    
## Printing the Diffusive Matrix
def printK(K,params):
    dim = params["dim"]
    ksmall = []
    tmp = []
    print("The diffusive matrix is:\n")
    for k in range (0,dim):
        ksmall = K[k]
        for j in range(0,dim):
      	    tmp = ksmall[j]
      	    #print(tmp)
      	    for i in range(0,dim+2):
                for l in range(0,dim+2):  
                    print("K[",k,",",j,",",i,",",l,"]=",tmp[i,l],"\n")
       
    return 0






