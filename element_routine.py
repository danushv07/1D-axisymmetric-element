# -*- coding: utf-8 -*-
"""
This is the element routine file wherein Jacobian, strain disp. mat,
elemental stiffness, elemental F_int are computed.

The material routine is called in F_e_int() function, which returns the stress,
overstress and material tangent stiffness .


@author: Danush
"""

import numpy as np
from material_routine import *

##############################################################################

# N- shape functions for the given problem
si = 0
N = np.array([(1/2)*(1-si),(1/2)*(1+si)]).reshape(2,1)


#computation of jacobian
def jacobian_of_elements(r_elements):
    """This function returns the Jacobian for the problem.
    
    Input---> r_elements, depending on the meshgenerator
    
    Output--> Jacobian"""
    
    J = (r_elements[0,1] - r_elements[0,0])/2
    return J
 

# B strain displacement matrix
def strain_disp_mat(r_elements):
    """This function returns the strain displacement matrix.
    
    Input----> r_elements,N-shape functions
    
    Output---> B-strain displacement matrix(2x2)"""
    jacobian = jacobian_of_elements(r_elements)
    B = np.array([-1/(2*jacobian),1/(2*jacobian),
                        N[0]/(N[0]*r_elements[0,0] + N[1]*r_elements[0,1]),
                        N[1]/(N[0]*r_elements[0,0] + N[1]*r_elements[0,1])]).reshape(2,2)
    
    return B

# The no. of Gauss points for the problem is 1 , so w is assigned the value 2
# gauss points = 1
w = 2


# Elemental stiffness matrix
def k_e_mat(B,r,C_t,jacobian,w):
    """This function returns the elemental stiffness matrix
    which is called inside F_e_int function. 
    
    Input----> B-strain disp.Mat,r-radius from mesh generator,
               C_t-Tangent stiffness Mat.,jacobian,w 
               
    Output---> K_elemental"""
    
    
    k = np.array([w *(B.transpose().dot(C_t.dot(B)))*(N.transpose().dot(r.transpose()))
                 * jacobian]).reshape(2,2)
    return k


#Elemental Assignment matrix
def a_e_mat(n_e):
   """This function returns the elemental assignment matrix
   for the required no. of elements.
    
   Input-----> n_e no. of elements
   
   Output----> Elemental Assignment matrix"""
       
   B = np.zeros((2,n_e+1))
   h =[]
   for i in range(0,n_e):
       B[0][i] = 1
       B[1][i+1] = 1
       h.append(B)
       B = np.zeros((2,n_e+1))
   return h

#Elemental F_int and K-Stiffness
def F_e_int(u_e,delta_u_e,r,pre_stress,delta_t,T,Q,zetta,E,mu,n_e):
    """This function computes the elemental F_int, elemental K,
    stress, strain.
    
    To obtain C_t and over stress --> ovestress function is called in the 
    material routine
    
    Input---> u_e,delta_e-nodal displacements,r- with r_elements, 
              pre_stress- previous over stress, delta_t -time increment
              T,Q,zetta,E,mu,n_e - input parameters file()
              
    Output--> strain,stress,over_stress,F_internal_elements"""
    
    B = strain_disp_mat(r)
    jacobian = jacobian_of_elements(r)
    
    #calculation of elemental strains   
    #calc of elemental delta strains
    strain_e = B.dot(u_e)
    strain_del_e = B.dot(delta_u_e)
    
    # function call to material routine
    stress_e,over_stress_e,C_t = ovestress(delta_t,T,Q,pre_stress,strain_del_e,strain_e,zetta,E,mu)
    
    #computation of  elemental F_int
    F_int_e = w *(B.transpose().dot(stress_e).dot((N.transpose().dot(r.transpose())))
                     *jacobian)
    
    # function call to compute elemental K-stiffness
    elemental_k = k_e_mat(B,r,C_t,jacobian,w)
    
    return strain_e,stress_e,over_stress_e,F_int_e,elemental_k    
    
 
        
   
    
        
    



