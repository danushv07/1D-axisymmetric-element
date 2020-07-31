# -*- coding: utf-8 -*-
"""
@author: Danush
   This file contains all the necessary input parameters for the infinite element
   problem.
   Any changes here reflects in the main.py file
"""


import numpy as np

def parameter_list():
    """This function returns all the input parameters
    n_e - no. of elements to be considered
    zetta - value for the Time intergration scheme"""
    # Input parameters
    E      = 200000 #Young's Modulus in MPa
    mu     = 0.20 #poison's ratio
    Q      = 100000 #MPa
    T      = 1 #sec
    a      = 50 #mm
    b      = 100 #mm
    p_max  = 140 #MPa
    t_l    = 2 #sec
    t_f    = 10 #sec
    n_e    = 10 # no. of elements
    zetta = 1/2 #Euler modified method
    del_t = 0.1 # delta time/time step in s
    return E,mu,Q,T,a,b,p_max,t_l,t_f,n_e,zetta,del_t

   