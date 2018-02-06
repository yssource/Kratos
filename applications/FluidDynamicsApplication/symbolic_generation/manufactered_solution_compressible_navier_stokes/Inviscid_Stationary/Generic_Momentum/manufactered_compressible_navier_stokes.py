## Manufactured Solutions for 2D Compressible Navier-Stokes problem
from KratosMultiphysics import *

import matplotlib.pyplot as plt
import numpy as np
from scipy import linspace, meshgrid, arange, empty, concatenate, newaxis, shape
from mpl_toolkits.mplot3d.axes3d import Axes3D

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

## Source Terms derived from manufactured solutions
#S_rho = Symbol("S_rho", real = True) NOT NEEDED IF THE DENSITY IS CONSTANT IN TIME AND MOMENTUM IS SOLENOIDAL
S_u = Symbol("S_u", real = True)
S_v = Symbol("S_v", real = True)
S_e = Symbol("S_e", real = True)
f_ext = DefineVector('f_ext',2)
r = Symbol('r', real=True)  

# Linear analytical expression of variables
R = 287.2
rho = 0.5*x+y+0.8     #ALWAYS POSITIVE (BETWEEN 0.8-3.5 ca)
u = 32
v = 8
m_u =  u*(1*x+2*y+5)         #DIVERGENCE FREE
m_v = v*(2*x-4*y+3)              #DIVERGENCE FREE
et = 100*x -500*y +36500
'''
Constant expressions of variables
rho = 2.7
m_u = 240
m_v = -12
et = 36500
'''
print("rho = ", rho)
print("m_U = ", m_u)
print("m_v = ", m_v)
print("et = ", et)

'''
print("CHECK DIVERGENCE: ", m_u, "      ", m_v)

for x in range(0,4):
    for y in range(0,2):

        #rho = 0.5*x+y+0.8
        u = 32
        v = 8
        m_u =  u*(1*x+2*y+5)
        m_v = v*(2*x-4*y+3)
        mod = np.sqrt(m_u**2+m_v**2)
        et = 100*x -500*y +36500
    
        T = (1.4-1)*(et-(m_u**2+m_v**2)/(2*rho))/(rho*R)
        Ma = mod/(rho*sqrt(1.4*R*T)) #  SUPERSONIC
        if rho>0:
            print("CHECK POSITIVE Density: ",rho)
        if et>0:
            print("CHECK POSITIVE Energy: ",et)
        if T<0:
            print("ERROR T")
        else:
            print("CHECK POSITIVE T: ", T)
            if Ma<1:
                print("ERROR Ma not Supersonic")
                print("SUBSONICA Mach: ", Ma)
            else:
                print("CHECK SUPERSONIC Mach: ",Ma)
        tmp = et/rho
        tmp -= (m_u**2+m_v**2)/(2*rho**2)   # MUST BE ALWAYS POSITIVE FOR THE SPEED OF SOUND
        
        if tmp<0:
            print("ERROR FOR x = ",x,"and y = ",y,"\n")
        else:
            print("CHECK SPEED OF SOUND: ",sqrt(tmp))
'''

# Definition of costants
mu = 0.0
print("mu", mu)
l = 0.0
print("lambda", l)
gamma = 1.4
cv = 718


p = (gamma-1)*(et-(m_u**2+m_v**2)/(2*rho))
'''
Tau_test using velocities
#tau_uxx = mu*(diff((m_u/rho),x)+diff((m_u/rho),x))-(2*mu)*(diff(m_u/rho,x)+diff(m_v/rho,y))/3
#tau_uyy = mu*(diff(m_v/rho,y)+diff(m_v/rho,y))-(2*mu)*(diff(m_u/rho,x)+diff(m_v/rho,y))/3
#tau_uxy = mu*(diff(m_u/rho,y)+diff(m_v/rho,x))
#print("tau_uxx = ", simplify(tau_uxx))
#print("\n\ntau_uyy = ", simplify(tau_uyy))
#print("\n\ntau_uxy = ", simplify(tau_uxy))
'''
tau_xx = (4*mu/(3*rho))*diff(m_u,x)-(2*mu/(3*rho))*diff(m_v,y)-(4*mu/(3*rho**2))*m_u*diff(rho,x)+(2*mu/(3*rho**2))*m_v*diff(rho,y)
tau_yy = (4*mu/(3*rho))*diff(m_v,y)-(2*mu/(3*rho))*diff(m_u,x)-(4*mu/(3*rho**2))*m_v*diff(rho,y)+(2*mu/(3*rho**2))*(m_u*diff(rho,x))
tau_xy = (mu/rho)*(diff(m_u,y)+diff(m_v,x))-(mu/rho**2)*(m_u*diff(rho,y)+m_v*diff(rho,x))
q_x = (l*et/(rho**2*cv))*diff(rho,x) \
    -(l/(rho**3*cv))*(m_u**2+m_v**2)*diff(rho,x) \
    +(l/(cv*rho**2))*(m_u*diff(m_u,x)+m_v*diff(m_v,x))-l*diff(et,x)/(rho*cv)
q_y = (l*et/(rho**2*cv))*diff(rho,y) \
    -(l/(rho**3*cv))*(m_u**2+m_v**2)*diff(rho,y) \
    +(l/(cv*rho**2))*(m_u*diff(m_u,y)+m_v*diff(m_v,y))-l*diff(et,y)/(rho*cv)

# Mass Conservation Equation
S_rho = diff(rho,t)+diff(m_u,x)+diff(m_v,y) #MUST BE ZERO
print("\nS_rho:\n",S_rho)

# Moment Equation: x direction
S_u = diff(m_u,t)+diff((m_u**2/rho+p-tau_xx),x)+diff((m_u*m_v/rho-tau_xy),y)
#print("\nS_u:\n",S_u)

# Moment Equation: y direction
S_v = diff(m_v,t)+diff((m_v**2/rho+p-tau_yy),y)+diff((m_u*m_v/rho-tau_xy),x)
#print("\nS_v:\n",S_v)

# Energy Equation
S_e = diff(et,t) \
    +diff( (  (et+p)*m_u/rho-m_u*tau_xx/rho-m_v*tau_xy/rho+q_x),x) \
    +diff( (  (et+p)*m_v/rho-m_u*tau_xy/rho-m_v*tau_yy/rho+q_y),y)
#print("\nS_e:\n",S_e)

f_ext[0] = S_u/rho
f_ext[1] = S_v/rho

print("\n\nf_ext0:", ccode(f_ext[0]))
print("\nf_ext1:\n",ccode(f_ext[1]))
#print("\n\nf_ext0:", factor(f_ext[0]))
#print("\nf_ext1:\n",factor(f_ext[1]))

r = (S_e-(f_ext[0]*m_u+f_ext[1]*m_v) )/rho
print("\nr:\n",ccode(r))
#print("\nr:\n",factor(r))


'''
## PLOT OF SOURCE TERMS
X = linspace(0,4.1,30)
Y = linspace(0,0.5,30)
x, y = meshgrid(X, Y)

# F_ext0
z = 

fig = plt.figure(figsize=(0,5))
ax = Axes3D(fig)
ax.plot_surface(x,y,z, rstride=1, cstride=1)

# F_ext1
z =  

fig = plt.figure(figsize=(0,5))
ax = Axes3D(fig)
ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap='hot')

# r
z = 

fig = plt.figure(figsize=(0,5))
ax = Axes3D(fig)
ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap='hot')

plt.show()
'''