# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 07:40:06 2020

@author: AJ
"""

import numpy as np
import matplotlib.pyplot as plt

ly = 3
lx = 4
f = np.zeros((5,ly,lx),dtype=float)
nodetype = np.zeros((ly,lx),dtype=float)
ex = np.array([0,1,0,-1,0])
ey = np.array([0,0,1,0,-1])

for j in range(0,lx):
    for i in range(0,ly):
        if nodetype[i,j] <=0:
            for k in range(1,3):
                nextX = j + ex[k] 
                nextY = i + ey[k] 
                if (nextX < lx and nextY < ly):
                    print(j,i,k,nextX,nextY)
                    if (nodetype[nextY, nextX] <= 0):
                        ftemp = f[k,nextY, nextX]
                        f[k,nextY, nextX] = f[k+2,i, j]
                        f[k+2,i, j] = ftemp