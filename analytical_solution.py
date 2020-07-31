# -*- coding: utf-8 -*-
"""
@author: Danush
This file contains the function to return the analytical solution for the 
considered visco-elastic problem.
"""

import numpy as np
from input_parameters import parameter_list

# function call to import parameters
E,mu,Q,T,a,b,p_max,t_l,t_f,n_e,C,zetta = parameter_list()
################################################################

def analytical(a,b,E,mu,p_max,n_e,r_ele):
    """This function returns the analytical solution of the given problem
    
        Output---> analytical u_r"""
    u_r_elastic = np.zeros((n_e+1,1))
    for i in range(0,n_e+1):
        z = (1 - 2*mu)*r_ele[i] + (b**2 / r_ele[i])
        y = (p_max/E) * (a**2/(b**2 - a**2))
        x = 1+mu
        u_r_elastic[i] = x*y*z
    return u_r_elastic

  