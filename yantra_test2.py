# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 07:56:08 2020

@author: AJ
"""

import numpy as np

corner = [0,0]
length=[15,15]
cord = []
dx = 3
d = 2
grid_type ='midway'
data = zip(corner,length)
for x,l in zip(corner,length):
    if grid_type == 'nodal':
        cord.append(np.linspace(x, l, round((l-x)/dx)+1))
        
    elif grid_type == 'midway':
        x= x-dx/2; l = l + dx/2
        cord.append(np.linspace(x, l, round((l-x)/dx)+1))
        
if d == 2:
    x, y = cord
    y = y[::-1]
    shape = (len(y),len(x))
    ny, nx = shape
    
nodetype = -1. * np.ones(shape, order='F') # F for fortran style