#This program makes the plot for the average velocity of the particles in a the Z direction for 
# a single case

# Author : Bruno Blais
# Last modified : December 3rd

#Python imports
import os
import math
import numpy
import matplotlib.pyplot as plt
import sys
dt = 5e-6 # time step of the simulation

#=====================
#   Main plot
#=====================

fname=sys.argv[1]

#INPUT
print ("R-> ", fname)
N, u,v,w, unorm, wStd, wMax, wMin = numpy.loadtxt(fname, unpack=True)

t = N * dt
#Single sphere unsteady solution using Euler finite difference scheme
Dt = 5e-6	# time step of Euler method
rhof = 1000.	# fluid density
rhos = 1500.	# solid density
g = -10.	# gravity
dp = 0.0050	# particle diameter
mu = 0.001	# viscosity of the fluid	
m = 4*numpy.pi/3 * (dp/2)**3 * (rhos)
V= 4*numpy.pi/3 * (dp/2)**3 

drho = rhos - rhof

niter=int(max(t)/Dt)
vt=numpy.zeros([niter])
T=numpy.zeros([niter])
for i in range(0,niter-1,1):
    Rep = (dp * abs(vt[i])*rhof/mu)+10**-(18)
    Cd = (0.63 + 4.8/math.sqrt(Rep))**2
    Cdf=  0.125 * Cd *math.pi* rhof * dp**2 
    Fd = - 0.125 * Cd *math.pi* rhof * dp**2 * abs(vt[i]) * vt[i]
    T[i+1]=T[i]+Dt
    vt[i+1] = vt[i] + Dt/m * (V*(drho)*g+Fd )

print ('Mass is : ' , m)
print ('m *g *drho' , (drho)*g*V)
print ('Cd ' , Cdf * vt[i])
print ('Archimedes force', rhof * V)
print ('vt ', vt[niter-2])
print ('alpha', Cdf *dt*abs(vt[i])/m)
print ('The final drag force is to be : ' , Fd)
print ('Final velocity difference is : ', vt[-1] - w[-1])

plt.figure(fname)
plt.errorbar(t,w,yerr=wStd,fmt='ro', label="Average velocity")
plt.plot(t,wMax,'go', label="Maximal velocity")
plt.plot(T,vt,'-', label="Stokes analytical solution")
plt.ylabel('Average settling velocity [m/s]')
plt.xlabel('time [s]')
plt.title('Average velocity of the particles')
plt.legend(loc=9)
plt.show()

