# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 08:19:48 2020

@author: AJ
"""
import numpy as np

class Domain(object):
    """
    Generic parent class for creating a simulation domain 
    """
    _signature = 'yantra._base.Domain'
    def __init__(self,corner,lengths,dx,grid_type='midway',nodetype=None):
        """
        Initialises Domain instance
        
        Parameters
        ----------
        corner: tuple or list
            corner of the domain for 2D its left bottom corner for 3D it is back left bottom corner
        lengths: tuple or list
            length of domain in x,y directions for 2D and x,y,z directions in 3D
        dx: float
            discretization in physical units
        grid_type: str, optional
             can be specified as `midway` or `nodal`. By default it is set to `midway`
             In `midway` nodes are shifted by half the discretization so that each
             node represents volume around it and boundary nodes are located outside
             the domain of interest thus locating the boundary of the domain  in between boundary node
             and first interior node. This is analogous to cell-centered finite volume meshing.
             In `nodal` the first node represents boundary node and start or end of the domain. Nodes
             are not assumed to represent a volume.             
        nodetype: ndarray, optional
            represents array that can be used to mark a given node with a number and all the
            nodes that are then associated with that number can assign same parameters. All the
            nodetype > 0 are considered inactive nodes in physics and nodetype <= 0 considered as active nodes.
        """
        #check input
        self.d =2
        try:
            assert(self.d == len(corner))
            assert(self.d == len(lengths))
        except AssertionError:
            ValueError("Corner or lengths argument should contain %s values"%self.d)
        self.corner = corner
        self.lengths = lengths
        try:
            assert(grid_type.lower()=='nodal' or grid_type.lower()=='midway')
        except AssertionError:
            ValueError("grid_type can either be nodal or midway, %s given"%grid_type)        
        self.grid_type = grid_type.lower()
        self.dx = dx
        self._process_input()
        if np.all(nodetype == None):
            self.nodetype = nodetype
        if nodetype is None:
            self.nodetype = -1. * np.ones(self.shape, order='F')
        else:
            self.nodetype = np.asfortranarray(nodetype)
            
    def _process_input(self):
         """
         processes input for domain
         """
         d,dx,corner,lengths,grid_type = self.d,self.dx,self.corner,self.lengths,self.grid_type
         cord = []
         for x,l in zip(corner,lengths):
             if grid_type == 'nodal':
                 cord.append(np.linspace(x, l, round((l-x)/dx)+1))
             elif grid_type == 'midway':
                 x= x-dx/2; l = l + dx/2
                 cord.append(np.linspace(x, l, round((l-x)/dx)+1))
         if d == 2:
             self.x, self.y = cord
             self.y = self.y[::-1]
             self.shape = (len(self.y),len(self.x))
             self.ny, self.nx = self.shape
         elif d == 3:
             self.x,self.y,self.z = cord
             self.y = self.y[::-1]
             self.z = self.z[::-1]
             self.shape = (len(self.z),len(self.y),len(self.x))
             self.nz,self.ny,self.nx = self.shape
             
    def _bounding_box(self,center,lengths):
        """get bounding box"""
        d,dx = self.d,self.dx
        if d == 2:
            cord = [self.x,self.y]
        elif d == 3:
            cord = [self.x,self.y,self.z]
        bb = []
        for i in range(d):
            mini=np.abs(cord[i]-((center[i]-lengths[i]/2.)-dx)).argmin()
            maxi=np.abs(cord[i]-((center[i]+lengths[i]/2.)+dx)).argmin() 
            if i==0:
                bb.append((mini,maxi+1))
            else:
                bb.append((maxi,mini+1))
        return bb
    
    def draw_rect(self, center, lengths, idx=1.):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the rectangle

        Parameters
        ----------
        center: tuple or list
            center of rectangle
        lengths: tuple or list
            length of rectangle in x and y direction
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the rectangle
        
        See also
        --------
        :func:`draw_circle`,`draw_multicoated_circle`
        """
        x,y = self.x, self.y
        lengths = (lengths[0] + self.dx ,lengths[1] + self.dx)
        bb = self._bounding_box(center,lengths)
        print(bb)
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        xx, yy = self.meshgrid(x[xmin:xmax],y[ymin:ymax])
        marker = 1. * ((xx >= (center[0] - lengths[0]/ 2.)) * 
                       (xx <= (center[0] + lengths[0]/ 2.)) *
                       (yy >= (center[1] - lengths[1]/ 2.)) * 
                       (yy <= (center[1] + lengths[1]/ 2.)))
        self.nodetype[ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[ymin:ymax,xmin:xmax] * (marker <= 0))
        return
    
    def meshgrid(self,*args):
         """
         creates co-ordinates in meshgrid format
         
         Parameters
         ----------
         x,y: list or 1D ndarray (for 2D), optional
         x,y,z: list or 1D ndarray (for 3D), optional
         
         Returns
         -------
         x,y: ndarray (for 2D)
         x,y,z: ndarray (for 3D)
         """
         if self.d == 2:
             if len(args) > 0 and len(args)<=2:
                 try:
                     x,y  = args
                 except IndexError:
                     ValueError("There should be two input lists")
             else:
                 x,y = self.x,self.y                    
             return np.meshgrid(x,y)
         elif self.d == 3:
             if len(args) > 0 and len(args)<=3:
                 try:
                     x,y,z  = args
                 except IndexError:
                     ValueError("There should be three input lists")
             else:
                 x,y,z = self.x,self.y,self.z
             y,z,x = np.meshgrid(y,z,x)
             return x,y,z