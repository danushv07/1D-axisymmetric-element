# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 18:49:40 2020

@author: Danush
"""
""" Generate list of position of nodes according to a geometric series
#    for assignement in "Nonlinear Finite Element Methods" 
#    in summer term 2020
#    lecturer in charge: Dr. Geralf Hütter"""

import matplotlib.pyplot as plt
import numpy as np
#####################################################################

def r_element(a,b,n_e):
    """This function returns the radial distance of each elements(nodal values)
    
    Input--> inner(a) and outer(b) radius, n_e->no. of elements
    
    Output-> node location"""
    
    #Input Parameters
    meshrefinementfactor = 2 #ratio of element sizes at outer and inner radius
    
    #ratio between element sizes of subsequent elements for a geometric series
    q = meshrefinementfactor**(1/(n_e-1))
    #size of first interval
    dr = (b-a)*(1-q)/(1-meshrefinementfactor*q)
    rnode = a
    rnodes = [a]
    #loop over all elements
    for i in range(0,n_e):
        rnode=rnode+dr
        rnodes.append(rnode)
        dr = dr*q
    return np.array(rnodes)
##########################################
#plotting    
#visualize location of nodes
#plt.plot(g,np.zeros((nelem+1)),'x')
#plt.xlabel('r')
#plt.legend('nodes')

