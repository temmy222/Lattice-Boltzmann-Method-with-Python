# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 20:20:55 2019

@author: tajayi3
"""
'''
#2D diffusion
import numpy as np
import matplotlib.pyplot as plt


n = 10
m = 10
xaxis = np.linspace(1,m,m)
yaxis = np.linspace(1,n,n)
q = 9
f = np.zeros((q,n,m),dtype=float)
rho = np.zeros((n,m),dtype=float)
x = np.zeros((n),dtype=float)
y = np.zeros((m),dtype=float)
dx = 1
dy=dx
dt = 1
x[0] =0
y[0] = 0
for i in range(1,n):
    x[i] = x[i-1]+dx
    y[i] = y[i-1]+dy
    
csq = dx*dx/(dt*dt)
tw = 1
alpha = 0.25  # diffusion coefficient
omega = 1/(2*alpha/(dt*csq)+0.5) # relaxation time/factor

mstep = 100 # total number of time steps
w0 = 4.0/9.0
w1= 1.0/9.0
w2 = 1.0/36.0

w = ([w0,w1,w1,w1,w1,w2,w2,w2,w2])
for j in range(0,m):
    for i in range(0,n):
        rho[i,j]=1

#collision
for j in range(0,m):
    for i in range(0,n):
        for k in range(0,q):
            feq = w[k]*rho[i,j]
            f[k,i,j]=omega*feq+(1.-omega)*f[k,i,j]
            
#streaming
for j in range(m-1,0,-1):
    for i in range(0,n-1):
        f[2,i,j] = f[2,i,j-1]
        f[6,i,j] = f[6,i+1,j-1]
for j in range(m-1,0,-1):
    for i in range(n-1,0,-1):
        f[1,i,j] = f[1,i-1,j]
        f[5,i,j] = f[5,i-1,j-1]
for j in range(0,m-1):
    for i in range(n-1,0,-1):
        f[4,i,j] = f[4,i,j+1]
        f[8,i,j] = f[8,i-1,j+1]
for j in range(0,m-1):
    for i in range(0,n-1):
        f[3,i,j] = f[3,i+1,j]
        f[7,i,j] = f[7,i+1,j+1]
        

#boundary condition
for j in range(0,m):
    f[1,0,j] = w[1]*tw + w[3]*tw -f[3,0,j]
    f[5,0,j] = w[5]*tw + w[7]*tw -f[7,0,j]
    f[8,0,j] = w[8]*tw + w[6]*tw -f[6,0,j]
    f[3,n-1,j] = -f[1,n-1,j]
    f[6,n-1,j] = -f[8,n-1,j]
    f[7,n-1,j] = -f[5,n-1,j]
    
for i in range(0,n):
    f[4,i,m-1] = -f[2,i,m-1]
    f[7,i,m-1] = -f[5,i,m-1]
    f[8,i,m-1] = -f[6,i,m-1]
    f[1,i,0] = f[1,i,1]
    f[1,i,0] = f[1,i,1]
    f[2,i,0] = f[2,i,1]
    f[3,i,0] = f[3,i,1]
    f[4,i,0] = f[4,i,1]
    f[5,i,0] = f[5,i,1]
    f[6,i,0] = f[6,i,1]
    f[7,i,0] = f[7,i,1]
    f[8,i,0] = f[8,i,1]
    
for j in range(0,m):
    for i in range(0,n):
        suma =0
        for k in range(0,q):
            suma = suma+f[k,i,j]
        rho[i,j] = suma

plt.plot(xaxis,rho[:,1])
plt.grid()
plt.figure()
plt.plot(yaxis,rho[1,:])
plt.grid()

'''
#2D Advection diffusion
import numpy as np
import matplotlib.pyplot as plt


n = 10
m = 10
xaxis = np.linspace(1,m,m)
yaxis = np.linspace(1,n,n)
q = 9
f = np.zeros((q,n,m),dtype=float)
rho = np.zeros((n,m),dtype=float)
x = np.zeros((n),dtype=float)
y = np.zeros((m),dtype=float)
feq = np.zeros(q)
dx = 1
dy=dx
u = 0.1
v = 0.4
dt = 1
x[0] =0
y[0] = 0
for i in range(1,n):
    x[i] = x[i-1]+dx
    y[i] = y[i-1]+dy
    
ck = dx/dt
csq = ck*ck
tw = 1
alpha = 1.0  # diffusion coefficient
omega = 1/(3*alpha/(dt*csq)+0.5) # relaxation time/factor

mstep = 50000 # total number of time steps
w0 = 4.0/9.0
w1= 1.0/9.0
w2 = 1.0/36.0

w = ([w0,w1,w1,w1,w1,w2,w2,w2,w2])
for j in range(0,m):
    for i in range(0,n):
        rho[i,j]=1

for rto in range(0,mstep):
    #collision
    for j in range(0,m):
        for i in range(0,n):
            feq[0] = w[0]*rho[i,j]
            feq[1] = w[1]*rho[i,j]
            feq[2] = w[2]*rho[i,j]
            feq[3] = w[3]*rho[i,j]
            feq[4] = w[4]*rho[i,j]
            feq[5] = w[5]*rho[i,j]
            feq[6] = w[6]*rho[i,j]
            feq[7] = w[7]*rho[i,j]
            feq[8] = w[8]*rho[i,j]
            for k in range(0,q):
                f[k,i,j]=omega*feq[k]+(1.-omega)*f[k,i,j]
            
    #streaming
    for j in range(m-1,0,-1):
        for i in range(0,n-1):
            f[2,i,j] = f[2,i,j-1]
            f[6,i,j] = f[6,i+1,j-1]
    for j in range(m-1,0,-1):
        for i in range(n-1,0,-1):
            f[1,i,j] = f[1,i-1,j]
            f[5,i,j] = f[5,i-1,j-1]
    for j in range(0,m-1):
        for i in range(n-1,0,-1):
            f[4,i,j] = f[4,i,j+1]
            f[8,i,j] = f[8,i-1,j+1]
    for j in range(0,m-1):
        for i in range(0,n-1):
            f[3,i,j] = f[3,i+1,j]
            f[7,i,j] = f[7,i+1,j+1]
        

#boundary condition
#left hand boundary condition given temperature as tw
    for j in range(0,m):
        f[1,0,j] = w[1]*tw + w[3]*tw -f[3,0,j]
        f[5,0,j] = w[5]*tw + w[7]*tw -f[7,0,j]
        f[8,0,j] = w[8]*tw + w[6]*tw -f[6,0,j]
    
#top boundary condition T=0
    for i in range(0,n):
        f[4,i,m-1] = -f[2,i,m-1]
        f[7,i,m-1] = -f[5,i,m-1]
        f[8,i,m-1] = -f[6,i,m-1]
        f[3,i,m-1] = -f[1,i,m-1]
        f[0,i,m-1] = 0.0
    
#right hand boundary condition T=0
    for j in range(0,m):
        f[4,i,n-1] = -f[2,i,n-1]
        f[7,i,n-1] = -f[5,i,n-1]
        f[8,i,n-1] = -f[6,i,n-1]
        f[3,i,n-1] = -f[1,i,n-1]
        f[0,i,n-1] = 0.0

## Adiabatic condition
#f[1,i,0] = f[1,i,1]
#f[2,i,0] = f[2,i,1]
#f[3,i,0] = f[3,i,1]
#f[4,i,0] = f[4,i,1]
#f[5,i,0] = f[5,i,1]
#f[6,i,0] = f[6,i,1]
#f[7,i,0] = f[7,i,1]
#f[8,i,0] = f[8,i,1]
    
    for j in range(0,m):
        for i in range(0,n):
            suma =0
            for k in range(0,q):
                suma = suma+f[k,i,j]
            rho[i,j] = suma

#plt.plot(xaxis,rho[:,1])
#plt.grid()
#plt.figure()
#plt.plot(yaxis,rho[1,:])
#plt.grid()
        
plt.figure()
X,Y = np.meshgrid(xaxis,yaxis)
cp = plt.contourf(X, Y, rho)