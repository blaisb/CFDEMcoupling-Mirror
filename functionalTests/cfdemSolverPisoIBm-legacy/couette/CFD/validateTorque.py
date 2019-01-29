# USAGE : python ./monitorTorque.py LOGFILE

# Author : Bruno Blais
# Last modified : 23-01-2014

#Python imports
#----------------
import os
import sys
import numpy
import time
import matplotlib.pyplot as plt
import re # Ouhh regular expressions :)
#----------------

#********************************
#   OPTIONS AND USER PARAMETERS
#********************************
omega = 1.*2.*numpy.pi
L=0.01
R=0.0238
k = 0.0064/0.0238
mu=1
#Analytical solution for the Torque in the couette geometry
# OK la solution analytique elle est bonne 
torque = -4.*numpy.pi*mu*omega*R*R*L*(k*k/(1.-k*k))
#torque= -4*numpy.pi*mu*omega*L*(1./((k*R)**(-2.) - R**(-2.))) 

#======================
#   MAIN
#======================

fname=sys.argv[1]

#Labeling
ax = plt.figure("Torque") #Create window
plt.ylabel('Torque [N-m]')
plt.xlabel('Time [s]')
plt.title('Dynamic evolution of the Torque')

[t,dragX,dragY,dragZ, momentX, momentY, momentZ] = numpy.loadtxt(fname, unpack=True)
plt.plot(t,momentZ,label='Numerical solution')
plt.plot(t,t/t*torque,'-',label='Analytical solution')
plt.legend()

print ("==============================================================================")
print ("                                ANALYSIS")
print ("Ratio of measured torque over analytical solution is : ", momentZ[-1]/torque)
print ("Ideal solution is 1")
print ("==============================================================================")
plt.show()

