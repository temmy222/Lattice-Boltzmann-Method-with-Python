# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:25:50 2019

@author: tajayi3
"""

import numpy as np
import matplotlib.pyplot as plt

"""
# 1D diffusion
m = 30
dt = 1
dx = 1
f0 = np.zeros(m)
f = np.zeros(m)
dx = 1
dt = 0.05

alpha = 0.25  # diffusion coefficient

mstep = 4000 # total number of time steps
twall = 1.0  # left hand wall temperature

for i in range(0,m):
    f0[i] = 0

f0[0] = twall
#f[0] = 1
f0[m-1] = f0[m-2]
#f[m-1] = f[m-2]
xaxis = np.linspace(1,m,m)
#
#for kk in range(1,mstep):
#    for i in range(1,m-1):
#        f[i] = f0[i] + dt*alpha*(f0[i+1]-2*f0[i-1]) / (dx*dx)
#    for i in range(1,m-1):
#        f0[i] = f[i] 
#    f0[m-1] = f0[m-2]
    

for kk in range(1,mstep):
    for i in range(1,m-1):
        f[i] = f0[i] + dt*alpha*(f0[i+1]-2*f0[i] + f0[i-1]) / (dx*dx)
    for i in range(1,m-1):
        f0[i] = f[i] 
    f0[m-1] = f0[m-2]
    f0[0]=twall
    
plt.plot(xaxis,f0)
plt.grid()
"""
# 1D Advection diffusion
u =0.1
m = 200
dt = 0.25
dx = 0.5
f0 = np.zeros(m) #old time step
f = np.zeros(m) #new time step


alpha = 0.25  # diffusion coefficient

mstep = 1600 # total number of time steps
twall = 1.0  # left hand wall temperature

for i in range(0,m):
    f0[i] = 0

f0[0] = twall
f0[m-1] = f0[m-2]
xaxis = np.linspace(1,m,m)


for kk in range(1,mstep):
    for i in range(1,m-1):
        adv = (dt*u*(f0[i+1] -f0[i]) ) /dx
        f[i] = (f0[i] + dt*alpha*(f0[i+1]-2*f0[i] + f0[i-1]) / (dx*dx)) - adv
    for i in range(1,m-1):
        f0[i] = f[i] 
    f0[m-1] = f0[m-2]
    f0[0]=twall
    
plt.plot(xaxis,f0)
plt.grid()