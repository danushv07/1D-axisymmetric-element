# -*- coding: utf-8 -*-
"""
This is the material routine for the given problem of visco-elastic material.

Time integration of constitutive equations are done to obtain and update the
stress, overstress, expressed in the zetta format, and C_t is also computed.
@author: Danush
"""


import numpy as np

########################################################################

# function to obtain the stress and over_stress
def ovestress(delta_t,T,Q,pre_stress,strain_del_e,strain_e,zetta,E,mu):
    """->The material tangent stiffness is computed. 
    
       ->The deviatoric strain is calculated to obtain the over_stress.
       
       ->Finally stress is computed.
       
       Input --> pre_stress-previous stress values, strain_e-elemental strain
             strain_del_e - elemental delta strain
             
       Output -> stress, over_stress, C_t"""
    C = (E/((1+mu)*(1-2*mu))) * np.array([(1-mu),mu,mu,(1-mu)]).reshape(2,2)
    C_t = C + (1/(1+((zetta*delta_t)/T))) * Q * np.array([2/3, -1/3, -1/3, 2/3]).reshape(2,2)
    
    # deviatoric strain
    alpha = 1 / (1+((zetta*delta_t)/T))
    betta = 1 - ((zetta*delta_t)/T)
    dev_strain = (np.array([2/3,-1/3,-1/3,2/3]).reshape(2,2)).dot(strain_del_e)
    
    #over_stress and stress
    stressov = alpha * ((betta * pre_stress) + (Q * dev_strain))
    stress = (C.dot(strain_e)) + stressov
    
    return stress,stressov,C_t


    







