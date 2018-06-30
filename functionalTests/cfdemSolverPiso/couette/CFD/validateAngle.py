
# Author : Bruno Blais

#Python imports
#==============================================
import os
import math
import numpy
import matplotlib.pyplot as plt
import sys
dt = 0.00001	# time step of the simulation
omega = 60./60. * 2. * 3.14159
R=0.0238
k = 0.0064/0.0238
omegaL=-omega
folder="results/"

#===============================================
#   Read the SRF variables
#===============================================

# read files
print ("R-> %s velocity")
t,x,y, vx, vy = numpy.loadtxt("particle.txt", unpack=True)

r = numpy.sqrt(x*x + y*y)
theta0= numpy.arctan(y/x)

#Analytical solution for theta velocity in eulerian frame of reference
u = omega *k* R * (-r/(R) + (R)/r) / (1/k - k)
uthx =- u * numpy.sin(theta0)
uthy = u * numpy.cos(theta0)


# Create the figures
#-----------------------------------------------

plt.figure("Speed comparison ")
plt.plot(t,vx,label = 'U" Simulation')
plt.plot(t,vy,label = 'V Simulation')
plt.plot(t,uthx,label = 'U Analytical')
plt.plot(t,uthy,label = 'V Analytical')
plt.legend(loc=3)
plt.ylabel('Variable')
plt.xlabel('timestep')

plt.show()


