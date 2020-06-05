# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 05:18:16 2020

@author: AJ
"""
from scipy.special import erfc
from numpy import exp
import numpy as np
import matplotlib.pyplot as plt

def analytical_soln(cb,x,ts,D, u):
    c = (cb/2)*(erfc((x - u*ts)/(4*D*ts)**0.5)+
         erfc((x+ u * ts)/(4*D *ts)**0.5)*exp(u*x/D))
    return c

cb = 1
ux=0.01
D=0.1
x = np.linspace(0,200,201)
time = 5000

cany = analytical_soln(cb,x,time,D,ux)
plt.plot(x,cany)

# class  Variable:
#     _args_name = ['type','default','dimension','isphyvar','skip_first_time','belongs_to','doc']
#     def __init__(self,*args):
#         for name,val in zip(self._args_name,args):
#             setattr(self,name,val)

# varos = {'c':Variable('scalar',0,{'N':1,'L':-3},True,False,'domain_params','concentration'),
#              'flux':Variable('vector',0,{'N':1,'L':-2,'T':-1},True,False,'domain_params','flux'),
#              'u':Variable('vector',0,{'L':1,'T':-1},True,False,'domain_params','velocity'),
#              'D':Variable('scalar',1./6.,{'L':2,'T':-1},True,False,'domain_params','diffusion coefficient'),
#              'Dref':Variable('parameter',1./6.,{'L':2,'T':-1},True,False,'solver_params','reference diffusion coefficient for setting timestep'),
#              'ss':Variable('scalar',0,{'N':1,'L':-3,'T':-1},True,False,'domain_params','source/sink term'),
#              'tau':Variable('scalar',1,{},False,False,'solver_params',
#                             'relaxation parameter/for TRT scheme its anti-symmetric component'),
#              'tauref':Variable('parameter',1,{},False,False,'solver_params',
#                             'refrence relaxation parameter/for TRT scheme its anti-symmetric component which is used to set timestep'),
#              'q':Variable('int',5,{},False,False,'solver_params','Number of lattice directions'),
#              'd': Variable('int',2,{},False,False,'solver_params','Dimension of the domain'),
#              'es2':Variable('parameter',1./3.,{},False,False,'solver_params','pseudo velocity of the sound'),
#              'magic_para':Variable('parameter',1./6.,{},False,False,'solver_params','Magic parameter for TRT LBM model'),
#              'collision_model': Variable('parameter','SRT',{},False,False,'solver_params','model for LB collision term'),
#              'interp': Variable('flag',1,{},False,False,'solver_params','interpolate to set correct boundary'),
#              'tfactbased': Variable('flag',0.,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/D'),
#              'tfact': Variable('parameter',1./6.,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/D'),
#              'lattice': Variable('parameter','D2Q5',{},False,False,'solver_params','lattice type for the model'),
#              'Dr': Variable('scalar',0.,{'L':2,'T':-1},True,False,'domain_params','remaing part of diffusion coefficient for diffusion velocity formulation'),
#              'u_adv': Variable('vector',0.,{'L':1,'T':-1},True,False,'domain_params','advective velocity specified in ade'),
#              'f': Variable('dist_func',0,{},False,True,'domain_params','LB distribution function'),
#              'time': Variable('parameter',0,{},False,True,'solver_params','simulation time in phyiscal units'),
#              'iters': Variable('parameter',0,{},False,True,'solver_params','number of iterations in simulations'),
#             }