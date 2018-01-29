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

# Analytical expression of variables
R = 287.2
#rho = 0.5*x+y+0.8     #ALWAYS POSITIVE (BETWEEN 0.8-3.5 ca)
rho = 2.7
u = 32
v = 8
m_u =  u*(1*x+2*y+5)         #DIVERGENCE FREE
m_v = v*(2*x-4*y+3)              #DIVERGENCE FREE
#m_u = 240
#m_v = -12
#et = 100*x -500*y +36500 #MUST BE POSITIVE
et = 36500
print("rho = ", rho)
print("m_U = ", m_u)
print("m_v = ", m_v)
print("et = ", et)


'''
print("CHECK DIVERGENCE: ", m_u, "      ", m_v)

for x in range(0,4):
    for y in range(0,2):

        #rho = 0.5*x+y+0.8
        rho = 2.7
        u = 32
        #print(u)
        v = 8
        m_u = 240
        m_v = -12
        #m_u =  u*(1*x+2*y+5)
        #m_v = v*(2*x-4*y+3)
        mod = np.sqrt(m_u**2+m_v**2)
        et = 100*x -500*y +36500
        #et = 36500
    
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
mu = 0.03
print("mu", mu)
l = 1.0
print("lambda", l)
gamma = 1.4
cv = 718


p = (gamma-1)*(et-(m_u**2+m_v**2)/(2*rho))
#tau_xx = mu*(4/3*diff((m_u/rho),x)-2/3*diff(m_v/rho,y)) #IT IS THE SAME
tau_xx = (4*mu/(3*rho))*diff(m_u,x)-(2*mu/(3*rho))*diff(m_v,y)-(4*mu/(3*rho**2))*m_u*diff(rho,x)+(2*mu/(3*rho**2))*m_v*diff(rho,y)
tau_yy = (4*mu/(3*rho))*diff(m_v,y)-(2*mu/(3*rho))*diff(m_u,x)-(4*mu/(3*rho**2))*m_v*diff(rho,y)+(2*mu/(3*rho**2))*(m_u*diff(rho,x))
tau_xy = (mu/rho)*(diff(m_u,y)+diff(m_v,x))-(mu/rho**2)*(m_u*diff(rho,y)+m_v*diff(rho,x))
q_x = (l*et/(rho**2*cv))*diff(rho,x)-(l/(rho**3*cv))*(m_u**2+m_v**2)*diff(rho,x)+(l/(rho**2*cv))*(m_u*diff(m_u,x)+m_v*diff(m_v,x))-l*diff(et,x)/(rho*cv)
q_y = (l*et/(rho**2*cv))*diff(rho,y)-(l/(rho**3*cv))*(m_u**2+m_v**2)*diff(rho,y)+(l/(rho**2*cv))*(m_u*diff(m_u,y)+m_v*diff(m_v,y))-l*diff(et,y)/(rho*cv)
#tau_xx = 0.0
#tau_xy = 0.0
##tau_yy = 0.0
#q_x = 0.0
#q_y = 0.0

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
S_e = diff(et,t)+diff(((et+p)*m_u/rho-m_u*tau_xx/rho-m_v*tau_xy/rho+q_x),x)+diff(((et+p)*m_v/rho-m_u*tau_xy/rho-m_v*tau_yy/rho+q_y),y)
#print("\nS_e:\n",S_e)

f_ext[0] = S_u/rho
f_ext[1] = S_v/rho
#print("\nf_ext0:\n",factor(f_ext[0]))
print("\n\nf_ext0:", ccode(factor(f_ext[0])))
print("\nf_ext1:\n",ccode(factor(f_ext[1])))



r = (S_e-(f_ext[0]*m_u+f_ext[1]*m_v))/rho

print("\nr:\n",ccode(factor(r)))
'''
## PLOT OF SOURCE TERMS
X = linspace(0,4.1,30)
Y = linspace(0,0.5,30)
x, y = meshgrid(X, Y)

# F_ext0
z = -0.111111111111111*(132.75*x**7 + 1858.5*x**6*y + 5403.6*x**6 + 10690.2*x**5*y**2 + 56986.56*x**5*y + 49872.24*x**5 + 32562.0*x**4*y**3 + 241963.2*x**4*y**2 + 435854.88*x**4*y + 210706.2*x**4 + 55908.0*x**3*y**4 + 520819.2*x**3*y**3 + 1480152.96*x**3*y**2 + 1484427.456*x**3*y + 484003.584*x**3 + 52344.0*x**2*y**5 + 579859.2*x**2*y**4 + 2410179.84*x**2*y**3 + 3830741.568*x**2*y**2 + 2581992.3456*x**2*y + 627573.10464*x**2 + 22608.0*x*y**6 + 290995.2*x*y**5 + 1836460.8*x*y**4 + 4252435.2*x*y**3 + 4504826.88*x*y**2 + 2252610.10944*x*y + 433458.118656*x + 2246.39999999999*y**7 + 35435.52*y**6 + 495659.52*y**5 + 1686024.576*y**4 + 2553713.0496*y**3 + 1990095.962112*y**2 + 784439.0240256*y + 124394.8253184)/((0.5*x + 1.0*y + 0.8)**6*(1.0*x + 2.0*y + 1.6)**2)
fig = plt.figure(figsize=(0,5))
ax = Axes3D(fig)
ax.plot_surface(x,y,z, rstride=1, cstride=1)

# F_ext1
z = -0.111111111111111*(689.85*x**7 + 5971.5*x**6*y + 7634.16*x**6 + 17397.0*x**5*y**2 + 52488.0*x**5*y + 33500.16*x**5 + 8837.99999999999*x**4*y**3 + 96321.6*x**4*y**2 + 163447.2*x**4*y + 73767.492*x**4 - 56052.0000000001*x**3*y**4 - 107481.6*x**3*y**3 + 61943.0400000002*x**3*y**2 + 191723.04*x**3*y + 81663.2064*x**3 - 126244.8*x**2*y**5 - 589766.4*x**2*y**4 - 871061.759999999*x**2*y**3 - 469086.624*x**2*y**2 - 27392.2560000003*x**2*y + 32518.434816*x**2 - 103824.0*x*y**6 - 720276.48*x*y**5 - 1677265.92*x*y**4 - 1816463.232*x*y**3 - 968731.5456*x*y**2 - 226401.73056*x*y - 12957.6075264*x - 29664.0*y**7 - 291456.0*y**6 - 917890.56*y**5 - 1403075.52*y**4 - 1174588.416*y**3 - 544222.49472*y**2 - 127797.16608*y - 11216.8009728)/((0.5*x + 1.0*y + 0.8)**6*(1.0*x + 2.0*y + 1.6)**2)
fig = plt.figure(figsize=(0,5))
ax = Axes3D(fig)
ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap='hot')

# r
z = 0.111111111111111*(1876765122.0*x**10 + 23089082702.4*x**9*y - 63938302295.4*x**9 + 128351247537.6*x**8*y**2 - 1172895707099.52*x**8*y - 1228800908084.49*x**8 + 484876304409.6*x**7*y**3 - 8876629673389.44*x**7*y**2 - 18740620748962.4*x**7*y - 8804661442182.96*x**7 + 1622861303616.0*x**6*y**4 - 36385870333478.4*x**6*y**3 - 120723116005349.0*x**6*y**2 - 116351930845578.0*x**6*y - 35114886793212.3*x**6 + 4832728185600.0*x**5*y**5 - 87551210790086.4*x**5*y**4 - 427420063774781.0*x**5*y**3 - 643138151885957.0*x**5*y**2 - 397598988705812.0*x**5*y - 87865282098847.8*x**5 + 10839029848012.8*x**4*y**6 - 121972239344794.0*x**4*y**5 - 901930470081722.0*x**4*y**4 - 1.92210955940701e+15*x**4*y**3 - 1.84085360786283e+15*x**4*y**2 - 830959475264704.0*x**4*y - 144119895079452.0*x**4 + 16235828722483.2*x**3*y**7 - 82652913395619.8*x**3*y**6 - 1.14160337838682e+15*x**3*y**5 - 3.33704892748786e+15*x**3*y**4 - 4.45171310153756e+15*x**3*y**3 - 3.09531007610314e+15*x**3*y**2 - 1.09377456693309e+15*x**3*y - 155451801761574.0*x**3 + 15015397825843.2*x**2*y**8 + 2034550356664.25*x**2*y**7 - 815032165073662.0*x**2*y**6 - 3.33498537624359e+15*x**2*y**5 - 5.91013828661146e+15*x**2*y**4 - 5.66809904212906e+15*x**2*y**3 - 3.07228737453247e+15*x**2*y**2 - 887782996821446.0*x**2*y - 106481250928374.0*x**2 + 7670601197568.0*x*y**9 + 37870000416921.6*x*y**8 - 270390367655784.0*x*y**7 - 1.74689052152098e+15*x*y**6 - 4.06145054743746e+15*x*y**5 - 5.09015589584126e+15*x*y**4 - 3.78104783025142e+15*x*y**3 - 1.67073533535682e+15*x*y**2 - 406758187091081.0*x*y - 42052786842382.1*x + 1631992264704.0*y**10 + 16307569867161.6*y**9 - 17980620652149.1*y**8 - 357371958291277.0*y**7 - 1.11820991874469e+15*y**6 - 1.78673365491928e+15*y**5 - 1.71722437656645e+15*y**4 - 1.03517743872642e+15*y**3 - 384480096619006.0*y**2 - 80581883473861.2*y - 7295005034606.26)/((0.5*x + 1.0*y + 0.8)**7*(1.0*x + 2.0*y + 1.6)**2*(359.0*x + 718.0*y + 574.4)**2)
fig = plt.figure(figsize=(0,5))
ax = Axes3D(fig)
ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap='hot')

plt.show()
'''