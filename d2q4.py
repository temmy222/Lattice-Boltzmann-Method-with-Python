# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 18:46:38 2019

@author: tajayi3
"""


# 2D Diffusion
import numpy as np
import matplotlib.pyplot as plt
n = 10 # x direction
m = 10 # y direction
f1 = np.zeros((n,m),dtype=float)
f2 = np.zeros((n,m),dtype=float)
f3 = np.zeros((n,m),dtype=float)
f4 = np.zeros((n,m),dtype=float)
rho = np.zeros((n,m),dtype=float)
feq = np.zeros((n),dtype=float)
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
alpha = 0.25  # diffusion coefficient
omega = 1/(2*alpha/(dt*csq)+0.5) # relaxation time/factor

mstep = 100 # total number of time steps
xaxis = np.linspace(1,m,m)
yaxis = np.linspace(1,n,n)

for j in range(0,m):
    for i in range(0,n):
        rho[i,j]=1
for j in range(0,m):
    for i in range(0,n):
        f1[i,j]=0.25*rho[i,j]
        f2[i,j]=0.25*rho[i,j]
        f3[i,j]=0.25*rho[i,j]
        f4[i,j]=0.25*rho[i,j]
        
for kk in range(1,mstep):
    for j in range(0,m):
        for i in range(0,n):
            feq=0.25*rho[i,j]
            f1[i,j]=omega*feq+(1-omega)*f1[i,j]
            f2[i,j]=omega*feq+(1-omega)*f2[i,j]
            f3[i,j]=omega*feq+(1-omega)*f3[i,j]
            f4[i,j]=omega*feq+(1-omega)*f4[i,j]
    for j in range(0,m):
        for i in range(0,n):
            f1[n-i-1,j] = f1[n-i-2,j]
            f2[n-i-2,j] = f2[n-i-1,j]
    for j in range(1,m):
        for i in range(0,n):
            f3[i,m-j-1] = f3[i,m-j-2]
            f4[i,m-j-2] = f4[i,m-j-1]
    for j in range(1,m):
        f1[0,j] = 0.5 - f2[0,j]
        f3[0,j] = 0.5 - f4[0,j]
        f1[n-1,j] = 0
        f2[n-1,j] = 0
        f3[n-1,j] = 0
        f4[n-1,j] = 0
    for i in range(1,n):
        f1[i,m-1] = 0
        f2[i,m-1] = 0
        f3[i,m-1] = 0
        f4[i,m-1] = 0
        f1[i,0] = f1[i,1]
        f2[i,0] = f2[i,1]
        f3[i,0] = f3[i,1]
        f4[i,0] = f4[i,1]
    for j in range(0,m):
        for i in range(0,n):
            rho[i,j] = f1[i,j]+f2[i,j]+f3[i,j]+f4[i,j]
        
plt.plot(xaxis,rho[:,1])
plt.grid()
plt.figure()
plt.plot(yaxis,rho[1,:])
plt.grid()


'''
# 2D Advection Diffusion

import numpy as np
import matplotlib.pyplot as plt
n = 200 # x direction
m = 10 # y direction
f1 = np.zeros((n,m),dtype=float)
f2 = np.zeros((n,m),dtype=float)
f3 = np.zeros((n,m),dtype=float)
f4 = np.zeros((n,m),dtype=float)
rho = np.ones((n,m),dtype=float)
x = np.zeros((n),dtype=float)
y = np.zeros((m),dtype=float)
dx = 1
u = 0.1
v =0.2
dy=dx
dt = 1
# x[0] =0
# y[0] = 0
# for i in range(1,n):
#     x[i] = x[i-1]+dx
#     y[i] = y[i-1]+dy
    
x = np.linspace(0,n,n+1)
y = np.linspace(0,m,m+1)
    
ck = dx/dt
csq = ck*ck
alpha = 0.25  # diffusion coefficient
omega = 1/(2*alpha/(dt*csq)+0.5) # relaxation time/factor

mstep = 1000 # total number of time steps
xaxis = np.linspace(1,m,m)
yaxis = np.linspace(1,n,n)

# for j in range(0,m):
#     for i in range(0,n):
#         rho[i,j]=1
for j in range(0,m):
    for i in range(0,n):
        f1[i,j]=0.25*rho[i,j]
        f2[i,j]=0.25*rho[i,j]
        f3[i,j]=0.25*rho[i,j]
        f4[i,j]=0.25*rho[i,j]
        

for kk in range(1,mstep):
    for j in range(0,m):
        for i in range(0,n):
            feq1=0.25*rho[i,j] * (1+2*(u/ck))
            feq2=0.25*rho[i,j] * (1-(2*u/ck))
            feq3=0.25*rho[i,j] * (1+2*(v/ck))
            feq4=0.25*rho[i,j] * (1-(2*v/ck))
            f1[i,j]=omega*feq1+(1-omega)*f1[i,j]
            f2[i,j]=omega*feq2+(1-omega)*f2[i,j]
            f3[i,j]=omega*feq3+(1-omega)*f3[i,j]
            f4[i,j]=omega*feq4+(1-omega)*f4[i,j]
    for j in range(0,m):
        for i in range(1,n):
            f1[n-i,j] = f1[n-i-1,j]
            f2[i-1,j] = f2[i,j]
    for i in range(0,n):
        for j in range(1,m):
            f3[i,m-j] = f3[i,m-j-1]
            f4[i,j-1] = f4[i,j]
    for j in range(1,m):
        f1[0,j] = 0.5 - f2[0,j]
        f3[0,j] = 0.5 - f4[0,j]
        f1[n-1,j] = 0
        f2[n-1,j] = 0
        f3[n-1,j] = 0
        f4[n-1,j] = 0
    for i in range(1,n):
        f1[i,m-1] = 0
        f2[i,m-1] = 0
        f3[i,m-1] = 0
        f4[i,m-1] = 0
        f1[i,0] = f1[i,1]
        f2[i,0] = f2[i,1]
        f3[i,0] = f3[i,1]
        f4[i,0] = f4[i,1]
    for j in range(0,m):
        for i in range(0,n):
            rho[i,j] = f1[i,j]+f2[i,j]+f3[i,j]+f4[i,j]

        
#plt.plot(xaxis,rho[:,1])
#plt.grid()
#plt.figure()
#plt.plot(yaxis,rho[1,:])
#plt.grid()

plt.figure()
X,Y = np.meshgrid(xaxis,yaxis)
cp = plt.contourf(X, Y, rho)
'''