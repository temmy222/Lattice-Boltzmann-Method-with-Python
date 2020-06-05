# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 00:24:36 2020

@author: AJ
"""
import numpy as np
import matplotlib.pyplot as plt

# D1Q2 for a constant temperature boundary condition
m=30
f0 = np.zeros((m),dtype=float)
f1 = np.zeros((m),dtype=float)
f2 = np.zeros((m),dtype=float)
rho = np.zeros((m),dtype=float)
feq = np.zeros((m),dtype=float)
x = np.zeros((m),dtype=float)
dt = 1
dx =1
thermal_cond = 0.34
flux = 30
for i in range(1,m):
    x[i] = x[i-1] +dx 

csq = (dx*dx)/(dt*dt)
alpha=0.25
omega= 1.0/(alpha/(dt*csq)+0.5)
mstep=200 # The total number of time steps
twall=1.0 # Left hand wall temperature
# Initial condition
for i in range(0,m):
    rho[i]=0.0 # Initial value of the domain temperature
    f1[i]=0.5*rho[i]
    f2[i]=0.5*rho[i]
 

# rho[0] = twall
# print(rho)

for kk in range(1,mstep):
    # main loop
    # collision process:
    for i in range(0,m):
        rho[i]=f1[i]+f2[i]
        feq[i]=0.5*rho[i]
        # since k1=k2=0.5, then feq1=feq2=feq
        f1[i]=(1-omega)*f1[i]+omega*feq[i]
        f2[i]=(1-omega)*f2[i]+omega*feq[i]
    # streaming process
    for i in range(1,m):
        f1[m-i]=f1[m-i-1] # f1 streaming
        f2[i-1]=f2[i] #f2 streaming
    # Boundary condition
    f1[0]=twall-f2[0]# constant temperature boundary condition, x=0
    # f1[m-1]=f1[m-2] # adiabatic boundary condition, x=L
    # f2[m-1]=f2[m-2]# adiabatic boundary condition, x=L
   
    
"""
Finite difference
"""
f_prev = np.zeros((m),dtype=float)
f = np.zeros((m),dtype=float)
f_prev[0]=twall
# f[0]=twall
# f_prev[m-1]=f_prev[m-2] # initial condition for old value of f at x=L
# f[m-1]=f[m-2] # initial condition for updated value of f at x=L



for kk in range(1,mstep):
    for i in range(1,m-1):
        f[i] = (f_prev[i] + dt*alpha*(f_prev[i+1]-2*f_prev[i] + f_prev[i-1]) / (dx*dx)) 
    for i in range(1,m):
        f_prev[i] = f[i] 
    # f_prev[m-1] = f_prev[m-2]
    
plt.plot(x,f_prev,marker='o',linestyle='--',label="FDM")
plt.plot(x,rho,marker='x',linestyle='--',label="LBM")
plt.grid()
plt.xlabel("x")
plt.ylabel("T")
plt.legend()