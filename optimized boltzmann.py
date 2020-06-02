# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 08:55:57 2020

@author: AJ
"""

import numpy as np
import matplotlib.pyplot as plt
c = np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1],[1,1],[-1,-1],[1,-1],[-1,1]])
ai = np.array([0,2,1,4,3,6,5,8,7])

na = 9  # number of lattice velocities
D=2 # Dimension of the simulation

w0 = 4.0/9.0
w1=1.0/9.0
w2 = 1.0/36.0

w = np.array([w0,w1,w1,w1,w1,w2,w2,w2,w2])

dt = 1; dx=1; S= dx/dt

c1 = 1.0
c2 = 3.0/(S**2)
c3 = 9.0/(2*S**4)
c4 = -3.0/(2.0*S**2)

# Initialize the relaxation time
nu_f = 0.1   # Viscosity

tau_f = nu_f * 3./(S*dt) +0.5

nt = 100 # Numer of time steps
nx = 101  # X-axis size
nz = 101 # Z-axis size

# Initialize arrays

f = np.zeros((na,nz,nx),dtype=float)
f_stream = np.zeros((na,nz,nx),dtype=float)
f_eq = np.zeros((na,nz,nx),dtype=float)
Delta_f = np.zeros((na,nz,nx),dtype=float)
rho = np.ones((nz,nx),dtype=float)
u = np.zeros((D,nz,nx),dtype=float)
Pi = np.zeros((D,nz,nx),dtype=float)
u2 = np.zeros((nz,nx),dtype=float)
cu = np.zeros((nz,nx),dtype=float)

# Initialize the density
rho_0 = 1.0
rho *=rho_0
rho[nz//2,3*nx//4] = 2*rho_0


for a in np.arange(na):
    f[a] = rho *w[a]
    
indexes = np.zeros((na,nx*nz),dtype=int)
for a in range(na):
    xArr = (np.arange(nx)-c[a][0]+nx)%nx
    zArr = (np.arange(nz)-c[a][0]+nz)%nz
    xInd, zInd = np.meshgrid(xArr,zArr)
    indTotal = zInd*nx +xInd
    indexes[a] = indTotal.reshape(nx*nz)
    
for t in np.arange(nt+1):
    #(1) periodic BC
    f[0:na,0:nz,0] = f[0:na,0:nz,-2]
    f[0:na,0:nz,-1] = f[0:na,0:nz,1]
    
    # (2) streaming step
    for a in np.arange(na):
        f_new = f[a].reshape(nx*nz)[indexes[a]]
        f_stream[a] = f_new.reshape(nz,nx)
    f=f_stream.copy()
    
    # (3) macroscopic properties: rho and u
    rho = np.sum(f,axis=0)
    Pi = np.einsum('azx,ad->dzx',f,c)
    u[0:D] = Pi[0:D]/rho
    
    # (4) Equiibirum distribution
    u2 = u[0]*u[0]+u[1]*u[1]
    for a in np.arange(na):
        cu = c[a][0]*u[0] +  c[a][1]*u[1]
        f_eq[a] = rho*w[a] * (c1+c2*cu+c3*cu**2+c4*u2)
    
    # (5) Collision term
    Delta_f = (f_eq-f)/tau_f
    f+=Delta_f
    
