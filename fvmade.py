# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:57:13 2020

@author: tajayi3
"""
import numpy as np
import matplotlib.pyplot as plt

u =0.1
m = 4
dt = 0.25
dx = 0.5
f0 = np.zeros(m) #old time step
f = np.zeros(m) #new time step
F1 = np.zeros(m) #Previous Residual
F2 = np.zeros(m) #Perturbed Residual
Jacobian = np.zeros(m) #Jacobian
jacMat = np.zeros((m,m))


alpha = 0.25  # diffusion coefficient

mstep = 3 # total number of time steps
twall = 1.0  # left hand wall temperature

for i in range(0,m):
    f0[i] = 0

#f0[0] = twall
rightb =0
#f0[m-1] = f0[m-2]
xaxis = np.linspace(1,m,m)
Cguess = f0
dP = np.ones(m)*1E-6

#for kk in range(1,mstep):
#    for i in range(1,m-1):
#        F1[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * f0[i+1] - 2*f0[i]+f0[i-1]/dx**2) + (u*(f0[i+1]-f0[i])/dx)
#        Cguess = f0+dP
#        print(F1)
#        F2[i] =  (Cguess[i] - f0[i])/dt - alpha * f0[i+1] - 2*f0[i]+f0[i-1]/dx**2 + u*(f0[i+1]-f0[i])/dx
#        #print(F2)
#        Jacobian[i] = F2[i]-F1[i]/dP[i]
#        jacMat[i-1,i-1] = Jacobian[i]
        
        
for kk in range(1,mstep):
    for i in range(0,m):
        if i==0:
            F1[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * f0[i+1] - 2*f0[i]+twall/dx**2) + (u*(f0[i+1]-f0[i])/dx)
        elif i==m-1:
            F1[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * rightb - 2*f0[i]+f0[i-1]/dx**2) + (u*(rightb-f0[i])/dx)
        else:
            F1[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * f0[i+1] - 2*f0[i]+f0[i-1]/dx**2) + (u*(f0[i+1]-f0[i])/dx)
    for i in range(0,m):
        Cguess = f0+dP
        if i==0:
            F2[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * f0[i+1] - 2*f0[i]+twall/dx**2) + (u*(f0[i+1]-f0[i])/dx)
        elif i==m-1:
            F2[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * rightb - 2*f0[i]+f0[i-1]/dx**2) + (u*(rightb-f0[i])/dx)
        else:
            F2[i] =  ((Cguess[i] - f0[i])/dt) - (alpha * f0[i+1] - 2*f0[i]+f0[i-1]/dx**2) + (u*(f0[i+1]-f0[i])/dx)
    for i in range(0,m):
        Jacobian[i] = F2[i]-F1[i]/dP[i]
        jacMat[i,i] = Jacobian[i]
    solu = np.linalg.solve(jacMat,F1)