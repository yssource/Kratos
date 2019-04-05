from sympy import *
from sympy.printing import cxxcode
import numpy as np

template_filename = "manufactured_template.cpp"
output_filename = "eca_flow_symbolic.cpp"

author_name = "Miguel Maso Sotomayor"

lower_case = "eca_flow"
cammel_case = ''.join(x.capitalize() for x in lower_case.split('_'))

# Constants
u = Symbol("mVelocity")
A = Symbol("mA")
B = Symbol("mB")
s = Symbol("mSigma")
w = Symbol("mOmega")
is_periodic = Symbol("mIsPeriodic")

# Coordinates
x = Symbol("x")
y = Symbol("y")
t = Symbol("rTime")

# Time dependent function
F = Piecewise((1.0 - cos(w*t), is_periodic), (1.0 - exp(-2.5*t), True))

# Auxiliary space dependent function
eta = s*y/x

# Velocity
u1 = u*erf(eta) + u*A*y*exp(-B*y)*sin((1 - 2*x)*pi) * F
u2 = u/s/sqrt(pi)*(1 - exp(-eta**2)) + 2*u*A*pi/B**2*cos((1 - 2*x)*pi)*(1 - exp(-B*y)*(B*y + 1)) * F

# Auxiliary functions for pressure
pa1 = x**3/3.0 - 3.0*x**2/4.0 + x/2.0 + 11.0/12.0
pa2 = y**2/2.0 + 7.0/8.0
pb1 = (4.0*x - 3.0)**4 - 2*(4.0*x - 3.0)**2 + 1.0
pb2 = 16.0*y**3 - 12.0*y**2 + 1

# Pressure
p = 50.0*log(pa1)*log(pa2) + 0.05*sin(pb1*pi/2.0)*sin(pb2*pi/2.0) * F

############################################################

## Velocity
# Time derivatives
du1_dt = diff(u1, t)
du2_dt = diff(u2, t)

# First derivatives
du1_dx1 = diff(u1, x)
du1_dx2 = diff(u1, y)
du2_dx1 = diff(u2, x)
du2_dx2 = diff(u2, y)

# Second derivatives
ddu1_dx11 = diff(du1_dx1, x)
ddu1_dx22 = diff(du1_dx2, y)
ddu2_dx11 = diff(du2_dx1, x)
ddu2_dx22 = diff(du2_dx2, y)

## Pressure
dp_dx1 = diff(p, x)
dp_dx2 = diff(p, y)

## Initialize the outstring to be filled with the template .cpp file
with open(template_filename) as template:
    outstring = template.read()

# Replace the names
outstring = outstring.replace("AuthorName", author_name)
outstring = outstring.replace("manufactured_template", lower_case)
outstring = outstring.replace("ManufacturedTemplate", cammel_case)

# Replace the computed derivatives
outstring = outstring.replace("// substitute U1", cxxcode(u1))
outstring = outstring.replace("// substitute U2", cxxcode(u2))
outstring = outstring.replace("// substitute DU1_DT", cxxcode(du1_dt))
outstring = outstring.replace("// substitute DU2_DT", cxxcode(du2_dt))
outstring = outstring.replace("// substitute DU1_DX1", cxxcode(du1_dx1))
outstring = outstring.replace("// substitute DU1_DX2", cxxcode(du1_dx2))
outstring = outstring.replace("// substitute DU2_DX1", cxxcode(du2_dx1))
outstring = outstring.replace("// substitute DU2_DX2", cxxcode(du2_dx2))
outstring = outstring.replace("// substitute DDU1_DX11", cxxcode(ddu1_dx11))
outstring = outstring.replace("// substitute DDU1_DX22", cxxcode(ddu1_dx22))
outstring = outstring.replace("// substitute DDU2_DX11", cxxcode(ddu1_dx11))
outstring = outstring.replace("// substitute DDU2_DX22", cxxcode(ddu1_dx22))
outstring = outstring.replace("// substitute P", cxxcode(p))
outstring = outstring.replace("// substitute DP_DX1", cxxcode(dp_dx1))
outstring = outstring.replace("// substitute DP_DX2", cxxcode(dp_dx2))

## Write the modified template
with open(output_filename, 'w') as output:
    output.write(outstring)
    output.close
