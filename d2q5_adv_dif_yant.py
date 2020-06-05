# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 05:13:35 2020

@author: AJ
"""

import numpy as np
import matplotlib.pyplot as plt

ly = 10
lx = 10
f = np.zeros((5,lx,ly),dtype=float)
c = np.zeros((lx,ly),dtype=float)
flux = np.zeros((2,lx,ly),dtype=float)
u = np.zeros((2,lx,ly),dtype=float)
tau  = np.zeros((lx,ly),dtype=float)
tau_a  = np.zeros((lx,ly),dtype=float)
ss  = np.zeros((lx,ly),dtype=float) # source term
nodetype = np.zeros((lx,ly),dtype=float)
ex = np.array([0,1,0,-1,0])
ey = np.array([0,0,1,0,-1])
MagicPara =0
w1=1/6
w2=2/6


# compute edf
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            #f1
            f[0,i, j] = w1 * c[i, j] 
            #f2
            ftemp = (1+3*u[0, i, j])
            f[1,i, j] = w2 * c[i, j] * ftemp
            #f3
            ftemp = (1+3*u[1, i, j])
            f[2,i, j] = w2 * c[i, j] * ftemp
            #f4
            ftemp = (1-3*u[0, i, j])
            f[3,i, j] = w2 * c[i, j] * ftemp
            #f5
            ftemp = (1-3*u[1, i, j])
            f[4,i, j] = w2 * c[i, j] * ftemp
        else:
            f[0,i, j] = 0
            f[1,i, j] = 0
            f[2,i, j] = 0
            f[3,i, j] = 0
            f[4,i, j] = 0


#collision srt (single relaxation time)
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            ssij = ss [i, j]
            omega = (1./tau[i, j])
            omega1 = 1-omega
            #compute equillibrium distribution for f1
            ftemp = w1 * c[i, j] 
            #collision step for f1
            f[0,i, j] = omega1 * f[0,i, j] + omega * ftemp + w1 * ssij
            #compute equillibrium distribution for f2
            ftemp =(1+3*u[0, i, j])
            ftemp = w2 * c[i, j] * ftemp
            #collision step for f2
            f[1,i, j] = omega1 * f[1,i, j] + omega * ftemp + w2 * ssij
            #compute equillibrium distribution for f3
            ftemp =(1+3*u[1, i, j])
            ftemp = w2 * c[i, j] * ftemp
            #collision step for f3
            f[2,i, j] = omega1 * f[2,i, j] + omega * ftemp + w2 * ssij
            #compute equillibrium distribution for f4
            ftemp =(1-3*u[0, i, j])
            ftemp = w2 * c[i, j] * ftemp
            #collision step for f4
            f[3,i, j] = omega1 * f[3,i, j] + omega * ftemp + w2 * ssij
            #compute equillibrium distribution for f5
            ftemp =(1-3*u[1, i, j])
            ftemp = w2 * c[i, j] * ftemp
            #collision step for f3
            f[4,i, j] = omega1 * f[4,i, j] + omega * ftemp + w2 * ssij
            # !f2 becomes f4 and f4 becomes f2
            #f3 becomes f5 and f5 becomes f3
            for k in range(1,3):
                ftemp = f[k+2,i, j]
                f[k+2,i, j] =f[k,i, j]
                f[k,i, j] = ftemp


                
#collision trt (two relaxation time)
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            ssij = ss [i, j]
            tau_s = 0.5 + (MagicPara/(tau_a[i, j]-0.5))
            omega_a = (1/tau[i, j])
            omega_s = 1 / tau_s
            #compute equilibirum distribution functions
            #f1
            f1eq = w1 * c[i, j] 
            #f2
            f2eq =(1+3*u[0, i, j])
            f2eq = w2 * c[i, j] * f2eq
            #f3
            f3eq =(1+3*u[1, i, j])
            f3eq = w2 * c[i, j] * f3eq
            #f4
            f4eq =(1-3*u[0, i, j])
            f4eq = w2 * c[i, j] * f4eq
            #f5
            f5eq =(1-3*u[1, i, j])
            f5eq = w2 * c[i, j] * f5eq
            # collision step
            #f1
            fseq = (f1eq+f1eq) * 0.5
            fs = (f[0,i, j]+f[0,i, j]) * 0.5
            f[0,i, j] = f[0,i, j] + omega_s * (fseq-fs) + w1 * ssij




# compute macro properties
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            c[i,j] = f[0,i,j]+f[1,i,j]+f[2,i,j]+f[3,i,j]+f[4,i,j]
            flux [0, i, j] = (1/(2*tau[i, j])) * c [i, j] * u [1,i, j] + (1-(1/(2*tau[i, j]))) * (f[1,i, j]-f[3,i, j])
            flux [1, i, j] = (1/(2*tau[i, j])) * c [i, j] * u [1,i, j] + (1-(1/(2*tau[i, j]))) * (f[2,i, j]-f[4,i, j])