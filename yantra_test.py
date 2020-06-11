# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 08:22:05 2020

@author: AJ
"""

from domain_yantra import Domain
domain = Domain((0,0),(15,15),3, grid_type = 'midway')
domain.draw_rect((7.5,7.5),(3,0.75),idx=1)
print(domain.nodetype)
pqty = 1.*(domain.nodetype>0)
domain_params={}
domain_params['D']=1e-9
domain_params['database']='cemdata07.dat'
domain_params['phrqc_input_file']='benchmark10.phrq'
domain_params['solution_labels']=100001
domain_params['eq_names']=['portlandite']
domain_params['solid_phases']={'portlandite':{'type':'non_diffusive','mvol':1,'c':pqty}}
domain_params['voxel_vol']=1
#solver parameters
solver_params={}
solver_params['collision_model']='srt'
solver_params['phrqc_flags']={}
solver_params['phrqc_flags']['only_interface']=True