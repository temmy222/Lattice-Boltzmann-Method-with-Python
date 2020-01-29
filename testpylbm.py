# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:27:28 2019

@author: tajayi3
"""

import numpy as np
import pylab as plt

def mesh(N):
    xmin, xmax = 0., 1.
    dx = 1./N
    x = np.linspace(xmin-.5*dx, xmax+.5*dx, N+2)
    return x

#x = mesh(10)
#plt.plot(x, 0.*x, 'sk')
#plt.show()

def equilibrium(m0, c):
    return c*m0

def initialize(mesh, c, la):
    m0 = np.zeros(mesh.shape)
    m0[np.logical_and(mesh<0.5, mesh>0.25)] = 1.
    print (m0)
    print (m0.shape)
    m1 = equilibrium(m0, c)
    print (m1)
    f0, f1 = np.empty(m0.shape), np.empty(m0.shape)
    print (f0)
    print (f1)
    m2f(f0, f1, m0, m1, la)
    return f0, f1, m0, m1

def f2m(f0, f1, m0, m1, la):
    m0[:] = f0 + f1
    m1[:] = la*(f1 - f0)

def m2f(f0, f1, m0, m1, la):
    f0[:] = 0.5*(m0-m1/la)
    f1[:] = 0.5*(m0+m1/la)

def relaxation(m0, m1, c, s):
    m1[:] = (1-s)*m1 + s*equilibrium(m0, c)

def transport(f0, f1):
    #periodical boundary conditions
    f0[-1] = f0[1]
    f1[0] = f1[-2]
    #transport
    f0[1:-1] = f0[2:]
    f1[1:-1] = f1[:-2]
    
# parameters
c = .5  # velocity for the transport equation
Tf = 2000. # final time
N = 12 # number of points in space
la = 1. # scheme velocity
s = 1.8 # relaxation parameter
# initialization
x = mesh(N)
f0, f1, m0, m1 = initialize(x, c, la)
t = 0
dt = (x[1]-x[0])/la
plt.figure(1)
plt.clf()
plt.plot(x[1:-1], m0[1:-1], 'k', label='init')
while t<Tf:
    t += dt
    relaxation(m0, m1, c, s)
    m2f(f0, f1, m0, m1, la)
    transport(f0, f1)
    f2m(f0, f1, m0, m1, la)
plt.plot(x[1:-1], m0[1:-1], 'r', label='final')
plt.legend()
plt.title('Advection')
plt.show()