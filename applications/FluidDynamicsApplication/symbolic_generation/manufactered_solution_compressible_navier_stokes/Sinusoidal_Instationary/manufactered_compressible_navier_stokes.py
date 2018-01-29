## Manufactured Solutions for 2D Compressible Navier-Stokes problem
from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint
 
t = Symbol("t", real = True)
x = Symbol("x", real = True)
y = Symbol("y", real = True)
gamma =  Symbol("gamma", real = True)       # NB IN THE FORMULATION IS 'y'!!!
mu =  Symbol("mu", real = True)             
l = Symbol("lambda", real = True)            
cv =  Symbol("cv", real = True)       

## Constant Variables for manufactured solutions
rho_0 = 0.5                 #Symbol("rho_0", real = True)
u_0 = 20.0                   #Symbol("u_0", real = True)
v_0 = 10.0                   #Symbol("v_0", real = True)
et_0 = 0.5                  #Symbol("et_0", real = True)
omg = 0.001                  #Symbol("omg", real = True)        # Omega
eps = 0.5                   #Symbol("eps", real = True)        # Epsilon
c = 1300                     #Inserted by me only for energy
mu = 0.3
l = 1.0
gamma = 1.4
R = 1.0
cp = 3.5
cv = 1.4
nu = 0.4                    #mu/rho_0

##  WITH THIS CHOICE MACH IS VERY SMALL, REMEMBER!!!

## Source Terms derived from manufactured solutions
S_rho = Symbol("S_rho", real = True)
S_u = Symbol("S_u", real = True)
S_v = Symbol("S_v", real = True)
S_e = Symbol("S_e", real = True)
f_ext = DefineVector('f_ext',2)
r = Symbol('r', real=True)  

# Analytical expression of variables
rho = rho_0*(sin(x**2+y**2+omg*t)+1.5)
u = u_0*(sin(x**2+y**2+omg*t)+eps)
v = v_0*(cos(x**2+y**2+omg*t)+eps)
m_u = rho*u
m_v = rho*v
#et = et_0*(cos(x**2+y**2+omg*t)+1.5)       # In the Salari Knupp test et = etot/rho of implemented formulation
et = et_0*(cos(x**2+y**2+omg*t)+c)       # In the Salari Knupp test et = etot/rho of implemented formulation
print("rho = ", rho)
print("m_u = ",m_u)
print("m_v = ",m_v)
print("et = ", et)


'''
for x in range(-1,1):
    for y in range(-1,1):
        for t in range(0,2):
            print("x = ",x," y = ",y," t = ",t)
            rho = rho_0*(sin(x**2+y**2+omg*t)+1.5)
            u = u_0*(sin(x**2+y**2+omg*t)+eps)
            v = v_0*(cos(x**2+y**2+omg*t)+eps)
            m_u = rho*u
            m_v = rho*v
            et = et_0*(cos(x**2+y**2+omg*t)+c)

            T = (gamma-1)*(et-(m_u**2+m_v**2)/(2*rho))/(rho*R)
            Ma = u/sqrt(gamma*R*T) #  SUPERSONIC
            if rho<0:
                print("CHECK POSITIVE Density: ", rho)
            if et<0:
                print("CHECK POSITIVE Energy: ",et)
            print("T: ", T)
            #if T<0:
                #print("CHECK POSITIVE T: ", T)
            if Ma<1:
                print("CHECK SUPERSONIC Mach: ",Ma)
            tmp = et/rho
            tmp -= (m_u**2+m_v**2)/(2*rho**2)   # MUST BE ALWAYS POSITIVE FOR THE SPEED OF SOUND
            if tmp<0:
                print("CHECK SPEED OF SOUND: ",sqrt(tmp))
                
'''


p = (gamma-1)*(et-(m_u**2+m_v**2)/(2*rho))
#tau_xx = mu*(4/3*diff((m_u/rho),x)-2/3*diff(m_v/rho,y)) #IT IS THE SAME
tau_xx = 4*mu/(3*rho)*diff(m_u,x)-2*mu/(3*rho)*diff(m_v,y)-4*mu/(3*rho**2)*m_u*diff(rho,x)+2*mu/(3*rho**2)*m_v*diff(rho,y)
tau_yy = 4*mu/(3*rho)*diff(m_v,y)-2*mu/(3*rho)*diff(m_u,x)-4*mu/(3*rho**2)*m_v*diff(rho,y)+2*mu/(3*rho**2)*(m_u*diff(rho,x))
tau_xy = mu/rho*(diff(m_u,y)+diff(m_v,x))-mu/rho**2*(m_u*diff(rho,y)+m_v*diff(rho,x))
q_x = l*et/(rho**2*cv)*diff(rho,x)-l/(rho**3*cv)*(m_u**2+m_v**2)*diff(rho,x)+l/(rho**2*cv)*(m_u*diff(m_u,x)+m_v*diff(m_v,x))-l*diff(et,x)/(rho*cv)
q_y = l*et/(rho**2*cv)*diff(rho,y)-l/(rho**3*cv)*(m_u**2+m_v**2)*diff(rho,y)+l/(rho**2*cv)*(m_u*diff(m_u,y)+m_v*diff(m_v,y))-l*diff(et,y)/(rho*cv)

# Mass Conservation Equation
S_rho = diff(rho,t)+diff(m_u,x)+diff(m_v,y) 
#print("\nS_rho:\n",S_rho)

# Moment Equation: x direction
S_u = diff(m_u,t)+diff((m_u**2/rho+p-tau_xx),x)+diff((m_u*m_v/rho-tau_xy),y)
#print("\nS_u:\n",S_u)

# Moment Equation: y direction
S_v = diff(m_v,t)+diff((m_v**2/rho+p-tau_yy),y)+diff((m_u*m_v/rho-tau_xy),x)
#print("\nS_v:\n",S_v)

# Energy Equation
S_e = diff(et,t)+diff(((et+p)*m_u/rho-m_u*tau_xx/rho-m_v*tau_xy/rho+q_x),x)+diff(((et+p)*m_v/rho-m_u*tau_xy/rho-m_v*tau_yy/rho+q_y),y)
#print("\nS_e:\n",S_e)

f_ext[0] = S_u/rho
f_ext[1] = S_v/rho
print("\nf_ext0:\n",factor(f_ext[0]))
print("\nf_ext1:\n",factor(f_ext[1]))

r = (S_e-(f_ext[0]*m_u+f_ext[1]*m_v))/rho

print("\nr:\n",factor(r))
