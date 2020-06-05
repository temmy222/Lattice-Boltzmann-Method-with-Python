# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:25:50 2019

@author: tajayi3
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from numpy import exp


#D1Q2 - Diffusion
#----------------------
'''
m = 30
dt = 1
dx = 1
f0 = np.zeros(m)
f1 = np.zeros(m)
f2 = np.zeros(m)
rho = np.zeros(m) # dependent variable
feq = np.zeros(m)
x = np.zeros(m)
x[0] = 0

for i in range(1,m):
    x[i] = x[i-1] + dx
    

csq = dx*dx/(dt*dt)
alpha = 0.25  # diffusion coefficient
omega = 1/(alpha/(dt*csq)+0.5) # relaxation time/factor

mstep = 200 # total number of time steps
twall = 1.0  # left hand wall temperature
f1[0] = twall
xaxis = np.linspace(1,m,m)

#initial condition
for i in range(0,m):
    rho[i] = 0   # Initial value of the domain temperature
    f1[i] = 0.5*rho[i] 
    f2[i] = 0.5*rho[i]

f1first = []
f1second = []
f2first = []
f2second = []
listt =[]
for kk in range(0,mstep):
    for i in range(0,m): # collision
        rho[i] = f1[i] +f2[i]
        feq[i] = 0.5 * rho[i] #0.5 is weighting factor
        f1[i] = (1-omega)*f1[i] + omega*feq[i]
        f2[i] = (1-omega)*f2[i] + omega*feq[i]
    for i in range(1,m): # streaming
        f1first.append(m-i)
        f1second.append(m-i-1)
        f2first.append(i-1)
        f2second.append(i)
        listt.append(i)
        f1[m-i] = f1[m-i-1] #  f1 streaming
        f2[i-1] = f2[i] # ! f2 streaming
    
        
    f1[0] = twall - f2[0] #constant temperature boundary condition, x=0
    f1[m-1] = f1[m-2] # adiabatic boundary condition, x=L
    f2[m-1] = f2[m-2] # adiabatic boundary condition, x=L

plt.plot(xaxis,rho)
plt.grid()
'''

#D1Q2 - Advection -  Diffusion
#----------------------

u =0.1
m = 100
dt = 1
dx = 1
f0 = np.zeros(m)
f1 = np.zeros(m)
f2 = np.zeros(m)
rho = np.zeros(m) # dependent variable
feq = np.zeros(m)
x = np.zeros(m)
x[0] = 0

for i in range(1,m):
    x[i] = x[i-1] + dx
    

ck = dx/dt
csq = ck*ck
alpha = 0.25  # diffusion coefficient
omega = 1/(alpha/(dt*csq)+0.5) # relaxation time/factor

mstep = 400 # total number of time steps
twall = 1.0  # left hand wall temperature
f1[0] = twall
xaxis = np.linspace(1,m,m)

#initial condition
for i in range(0,m):
    rho[i] = 0   # Initial value of the domain temperature
    f1[i] = 0.5*rho[i] 
    f2[i] = 0.5*rho[i]

f1first = []
f1second = []
f2first = []
f2second = []
listt =[]
for kk in range(0,mstep):
    for i in range(0,m): # collision
        rho[i] = f1[i] +f2[i]
        feq1 = 0.5 * rho[i] * (1+(u/ck)) #0.5 is weighting factor
        feq2 = 0.5 * rho[i] * (1-(u/ck)) #0.5 is weighting factor
        f1[i] = (1-omega)*f1[i] + omega*feq1
        f2[i] = (1-omega)*f2[i] + omega*feq2
    for i in range(1,m): # streaming
        f1first.append(m-i)
        f1second.append(m-i-1)
        f2first.append(i-1)
        f2second.append(i)
        listt.append(i)
        f1[m-i] = f1[m-i-1] #  f1 streaming
        f2[i-1] = f2[i] # ! f2 streaming
    
        
    f1[0] = twall - f2[0] #constant temperature boundary condition, x=0
    # f1[m-1] = 0 - f2[m-1] # boundary condition at x=L=0
#    f1[m-1] = f1[m-2] # adiabatic boundary condition, x=L
#    f2[m-1] = f2[m-2] # adiabatic boundary condition, x=L

plt.plot(xaxis,rho)
plt.grid()

#analytical solution
def analytical_soln(cb,x,ts,D, u):
    c = (cb/2)*(erfc((x - u*ts)/(4*D*ts)**0.5)+
         erfc((x+ u * ts)/(4*D *ts)**0.5)*exp(u*x/D))
    return c

c = np.zeros(101)
cb = 1
ts = 400
D =0.25
u = 0.1
x = np.linspace(0,100,101)
for i in range(0,100):
    c[i] = (cb/2)*(erfc((x[i] - u*ts)/(4*D*ts)**0.5)+erfc((x[i]+ u * ts)/(4*D *ts)**0.5)*exp(u*x[i]/D))
    
#plt.figure()
#plt.plot(x,c,xaxis,rho)
#plt.grid()

#finite difference solution
u =0.1
m = 100
dt = 0.25
dx = 1
f0 = np.zeros(m) #old time step
f = np.zeros(m) #new time step



alpha = 0.25  # diffusion coefficient

mstep = 1600 # total number of time steps
twall = 1.0  # left hand wall temperature
Peclet = alpha/(u*m)

for i in range(0,m):
    f0[i] = 0

f0[0] = twall
#f0[m-1] = f0[m-2]
f0[m-1] = 0
xaxis2 = np.linspace(1,dx*m,m)


for kk in range(0,mstep):
    for i in range(1,m-1):
#        adv = (dt*u*(f0[i+1] -f0[i-1]) ) /2*dx
        adv = (dt*u*(f0[i+1] -f0[i]) ) /dx  #upwind advection
        f[i] = (f0[i] + dt*alpha*(f0[i+1]-2*f0[i] + f0[i-1]) / (dx*dx)) - adv
    for i in range(1,m-1):
        f0[i] = f[i] 
#    f0[m-1] = f0[m-2]
    f0[m-1] = 0
    f0[0]=twall


plt.figure()
plt.plot(x,c,xaxis,rho,xaxis2,f0)
labels =['LBM','Analytical','Finite Difference']
plt.legend(labels)
plt.grid()
"""
import numpy as np
import matplotlib.pyplot as plt
c = np.array([[0,0],[1,0],[-1,0],[0,1],[1,1],[-1,-1],[1,-1],[-1,1]])
ai = np.array([0,2,1,4,3,6,5,8,7])

na = 9 #number of lattice velocities
D = 2 # Dimension of the simulation

w0 = 4.0/9.0
w1= 1.0/9.0
w2 = 1.0/36.0

w = ([w0,w1,w1,w1,w1,w2,w2,w2,w2])
dt=1; dx=1; S= dx/dt

c1 = 1.0
c2 = 3.0/(S**2)
c3 = 9.0/(2.0*S**4)
c4 = -3.0/(2*S**2)

nu_f = 0.1 # visocity
tau_f = nu_f *3./(S*dt) +0.5

nt = 100
nx = 11
nz = 11

f = np.zeros((na,nz,nx),dtype=float)
f_stream = np.zeros((na,nz,nx),dtype=float)
f_eq = np.zeros((na,nz,nx),dtype=float)
Delta_f = np.zeros((na,nz,nx),dtype=float)
rho = np.ones((nz,nx),dtype=float)
u = np.zeros((D,nz,nx),dtype=float)
Pi = np.zeros((D,nz,nx),dtype=float)
u2 = np.zeros((nz,nx),dtype=float)
cu = np.zeros((nz,nx),dtype=float)

rho_0 = 1.0
rho *= rho_0
rho[nz//2,3*nx//4] = 2*rho_0
    
for a in np.arange(na):
    f[a] = rho * w[a]
    
indexes = np.zeros((na,nz*nx),dtype=int)

"""

"""
u=0.1
dt=1
dx=1
D = 0.25
cb = 1
cfinite = np.zeros(11)

cforward = np.zeros(11)
for i in range(0,400):
    for j in range(1,len(cfinite)-1):
        cforward[j] = cfinite[j] - ((u*dt/2*dx)*(cfinite[j+1]-cfinite[j])) + D*dt/dx**2 * (cfinite[j+1] - 2*cfinite[j] +cfinite[j-1])
    cfinite=cforward
    cfinite[0] = cb
    cfinite[len(cfinite)-1] = 0 
"""