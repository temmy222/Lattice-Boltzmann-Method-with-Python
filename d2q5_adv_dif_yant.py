# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 05:13:35 2020

@author: AJ
"""

import numpy as np
import matplotlib.pyplot as plt

ly = 10 #essentially x
lx = 10 #essentially y
f = np.zeros((5,ly,lx),dtype=float)
c = np.zeros((ly,lx),dtype=float)
flux = np.zeros((2,ly,lx),dtype=float)
u = np.zeros((2,ly,lx),dtype=float)
u_adv = np.zeros((2,ly,lx),dtype=float)
tau  = np.zeros((ly,lx),dtype=float)
Dr  = np.zeros((ly,lx),dtype=float)
tau_a  = np.zeros((ly,lx),dtype=float)
ss  = np.zeros((ly,lx),dtype=float) # source term
nodetype = np.zeros((ly,lx),dtype=float)
ex = np.array([0,1,0,-1,0])
ey = np.array([0,0,1,0,-1])
MagicPara =0
topbc = 'periodic' 
bottombc = 'periodic'
leftbc = 'periodic' 
rightbc = 'periodic'
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
            # #f2 becomes f4 and f4 becomes f2
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
            #f2
            fseq = (f2eq+f4eq) * 0.5
            faeq = (f2eq-f4eq) * 0.5
            fs = (f[1,i, j]+f[3,i, j]) * 0.5
            fa = (f[1,i, j]-f[3,i, j]) * 0.5
            f[1,i, j] = f[1,i, j] + omega_s * (fseq-fs) + omega_a * (faeq-fa) + w2 * ssij
            #f4
            faeq = -faeq
            fa = -fa
            f[3,i, j] = f[3,i, j] + omega_s * (fseq-fs) + omega_a * (faeq-fa) + w2 * ssij
            #f3
            fseq = (f3eq+f5eq) * 0.5
            faeq = (f3eq+f5eq) * 0.5
            fs = (f[2,i, j]+f[4,i, j]) * 0.5
            fa = (f[2,i, j]-f[4,i, j]) * 0.5
            f[2,i, j] = f[2,i, j] + omega_s * (fseq-fs) + omega_a * (faeq-fa) + w2 * ssij
            #f5
            faeq = -faeq
            fa = -fa
            f[4,i, j] = f[4,i, j] + omega_s * (fseq-fs) + omega_a * (faeq-fa) + w2 * ssij
            # #f2 becomes f4 and f4 becomes f2
            #f3 becomes f5 and f5 becomes f3
            for k in range(1,3):
                fa = f[k+2,i, j]
                f[k+2,i, j] =f[k,i, j]
                f[k,i, j] = fa

def ud(Dr, tau, c, ua, fi, fii):
    Dtau = 3 * Dr / tau
    Dtau = Dtau / (1+Dtau)
    rame = Dtau * ((((fi-fii)/((c+1e-30)*(c+1e-30)))*c)-ua)
    return rame

#collision diffuse_velocity (diffuse_velocity formulation)
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            ssij = ss [i, j]
            omega = (1./tau[i, j])
            omega1 = 1-omega
            #update veolcity field using diffuse velocity formulation
            u[0, i, j] = ud (Dr[i, j], tau[i,j], c[i, j], u_adv[0, i, j], f[1,i, j], f[3,i, j])
            u[0, i, j] = u[0, i, j] + u_adv[0, i, j]
            u[1, i, j] = ud (Dr[i, j], tau[i,j], c[i, j], u_adv[1, i, j], f[2,i, j], f[4,i, j])
            u[1, i, j] = u[1, i, j] + u_adv[1, i, j]
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
            # #f2 becomes f4 and f4 becomes f2
            #f3 becomes f5 and f5 becomes f3
            for k in range(1,3):
                ftemp = f[k+2,i, j]
                f[k+2,i, j] =f[k,i, j]
                f[k,i, j] = ftemp


#stream and bounce-back
#Uses two step swap algorithim porposed by J. Latt(2008)
#performs streaming along with bounce back condition
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            for k in range(1,3):
                nextX = j + ex[k] 
                nextY = i + ey[k] 
                if (nextX < lx and nextY < ly):
                    if (nodetype[nextY, nextX] <= 0):
                        ftemp = f[k,nextY, nextX]
                        f[k,nextY, nextX] = f[k+2,i, j]
                        f[k+2,i, j] = ftemp


#Imposes boundary conditions on straight walls
#Accepted inputs
#topbc    = open,c,flux,periodic,nothing
#bottombc = open,c,flux,periodic,nothing
#leftbc   = open,c,flux,periodic,nothing
#rightbc  = open,c,flux,periodic,nothing
#location = midway/nodal
                        
if (topbc == 'periodic' or bottombc == 'periodic'):
    for j in range(0,lx):
        ftemp = f[4, j, 0]
        f[4, j, 0] = f[2, ly-1, j]
        f[2, ly-1, j] = ftemp
        
if (leftbc == 'periodic' or rightbc == 'periodic'):
    for i in range(0,ly):
        ftemp = f[3, j, 0]
        f[3, j, 0] = f[1, lx-1, j]
        f[1, lx-1, j] = ftemp
    

# compute macro properties
for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            c[i,j] = f[0,i,j]+f[1,i,j]+f[2,i,j]+f[3,i,j]+f[4,i,j]
            flux [0, i, j] = (1/(2*tau[i, j])) * c [i, j] * u [1,i, j] + (1-(1/(2*tau[i, j]))) * (f[1,i, j]-f[3,i, j])
            flux [1, i, j] = (1/(2*tau[i, j])) * c [i, j] * u [1,i, j] + (1-(1/(2*tau[i, j]))) * (f[2,i, j]-f[4,i, j])