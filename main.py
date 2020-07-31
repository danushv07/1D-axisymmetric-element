# -*- coding: utf-8 -*-
"""
@author: Danush
This file contains the function convergence_study() to solve the non linear problem.
-> The user is prompted to chose the case for linear or non linear behaviour.
-> Nodal displacements and stress, strain distributions are stored in .csv file.
-> Plots are generated after execution
"""

import numpy as np
import matplotlib.pyplot as plt
from input_parameters import *
from mesh_generator import *
from element_routine import *
from analytical_solution import *
import time
########################################################################
E,mu,Q,T,a,b,p_max,t_l,t_f,n_e,zetta,delta_t = parameter_list() #import all parameters
r_elements = r_element(a,b,n_e) # import nodal locations
u_r_elastic = analytical(a,b,E,mu,p_max,n_e,r_elements) #import analytical soln
################################################################################
print('NonLinear Finite Element Method Assignment 2020')
print('This is a program to solve visco-elastic material behaviour')
print('All the input parameters are specified in the input_parameters.py file')
print('Select either of the options')
lin = int(input('Linear elastic behaviour: 1 ;Nonlinear behaviour select: 0 :: '))

if lin == 1:
    Q = 0       
elif lin == 0:
    Q = Q
else:
    print('Choose the correct option!!! Run the program again')

def convergence_study():
    """This is the main function to compute displacements, stress, strain 
    
    Output---> nodal displacements, stress_distributions"""
    
    #initialization of variables
    time_seq = np.arange(0,10.1,delta_t) #from 0 - t_f in time steps
    F_ext = np.zeros((n_e+1,1),dtype=float) #initialize Global F_ext, dtype for linalg solve
    u = np.zeros((n_e+1,1)) # nodal displacemnets
    delta_u = np.zeros((n_e+1,1)) # delta displacements
    pre_stress = np.zeros((2,n_e)) # previous over stress
    u_rb_1 = [];u_rb_4 = [];u_rb_7 = []; u_rb_11= [] # list to obtain nodal disp. values for each time seq
    
    #time loop
    for i in time_seq:
        count = 0 #counter to check for 1 single iteration convergence
        if i<t_l:
            l_s = i * (1/t_l) #load scaling parameter time step*(1/t_l)
            F_ext[0] = p_max * a * l_s 
        else:
            F_ext[0] = p_max * a
        
        # Newton Raphson    
        while True:
            F_int = np.zeros((n_e+1,1),dtype=float) # initialize Global F_int
            K = np.zeros((n_e+1,n_e+1)) #initialize Global K-stiffness mat
            
            elemental_a = a_e_mat(n_e) # func. call to obtain Assignment mat
            stress_mat = [] # list to append stress values
            strain_mat = [] #list to append strain values
            for i in range(0,n_e):
                # elemental call
                u_e = np.array([u[i],u[i+1]]).reshape(2,1) 
                delta_u_e = np.array([delta_u[i],delta_u[i+1]]).reshape(2,1)
                r = np.array([r_elements[i],r_elements[i+1]]).reshape(1,2)
                
                # stre -  contains the pre. over stress values in 2x1 form
                stre = np.array([pre_stress[0,i],pre_stress[1,i]]).reshape(2,1)
                
                # element rountine call, then call the material routine internally
                strain,stress,over_stress,f_int,k_ele = F_e_int(u_e,delta_u_e,r,stre,delta_t,T,Q,zetta,E,mu,n_e)                
                
                #over stress assigned to stre and updated in pre. stress
                stre = over_stress
                pre_stress[0,i] = stre[0]
                pre_stress[1,i] = stre[1]
                stress_mat.append(stress)
                strain_mat.append(strain)
                # Global assigment of K and F_int
                K = K + (elemental_a[i].transpose().dot(k_ele.dot(elemental_a[i])))
                F_int = F_int +  elemental_a[i].transpose().dot(f_int)
            K = np.array(K, dtype='float')
            Residual = np.array((F_ext - F_int),dtype='float') # compute Residual
            delta_u = np.linalg.solve(K,Residual) #compute delta u
            
            u = u + delta_u #increment the u value
            count += 1
            
            # check for convergence
            if (np.linalg.norm(Residual,np.inf)<= 0.005*np.linalg.norm(F_int,np.inf) or np.linalg.norm(delta_u,np.inf) <= 0.005*np.linalg.norm(u,np.inf)):
                break
        # append only the last disp. values for plotting purpose
        u_rb_1.append(u[0]) ;u_rb_4.append(u[3]) ;u_rb_7.append(u[6]) ;u_rb_11.append(u[-1]) 
    
    return u,stress_mat,time_seq,strain_mat,u_rb_1,u_rb_4,u_rb_7,u_rb_11

# usage of time module and the main function call
start_time = time.time()   
convergence = convergence_study()
end_time = time.time()
print('total time taken: ',end_time-start_time)


#assignment of the values from function to respective variables
nodal_displacements = convergence[0]
stress              = convergence[1]
time_fol            = convergence[2]
strain              = convergence[3]
u_wide_1            = convergence[4]
u_wide_4            = convergence[5]
u_wide_7            = convergence[6]
u_wide_11           = convergence[7]


stress_rr = [] # list to get stress_rr values 
stress_phi = [] # list to get stress_phi_phi values 
strain_rr = [] # list to obtain strain_rr
strain_phi = [] # list to obtain strain_phi_phi
nodes = [0]
for i in range(0,n_e):
    nodes.append(i+1)
    stress_rr.append(stress[i][0])
    stress_phi.append(stress[i][1])
    strain_rr.append(strain[i][0])
    strain_phi.append(strain[i][1])
data = np.column_stack((stress_rr,stress_phi,strain_rr,strain_phi))  
data1 = np.column_stack((nodes,nodal_displacements)) 
    
np.savetxt('Stress_Strain distributions.csv',data) 
np.savetxt('nodal_displacements.csv',data1)


# Plotting
if lin == 1:
    plt.figure(1)
    plt.plot(r_elements,u_r_elastic,c='b',label='u_exact')
    plt.scatter(r_elements,nodal_displacements,c='r',label='u_numerical')
    plt.title('Convergence study')
    plt.xlabel('radius[mm]')
    plt.ylabel('$u_r$[mm]')
    plt.legend()
    plt.show()
    
elif lin == 0:
    plt.figure(1)
    plt.plot(r_elements,u_r_elastic,c='b',label='u_exact')
    plt.scatter(r_elements,nodal_displacements,c='r',label='u_numerical')
    plt.title('Convergence study')
    plt.xlabel('radius[mm]')
    plt.ylabel('$u_r$[mm]')
    plt.legend()
    plt.show()
    
    plt.figure(2)
    plt.plot(r_elements[:-1],stress_rr,c='r', marker='*', linestyle='--')
    plt.title('$\sigma_{rr}$ distribution')
    plt.xlabel('radius[mm]')
    plt.ylabel('$\sigma_{rr}$[MPa]')
    plt.show()
    
    
    plt.figure(3)
    plt.plot(r_elements[:-1],stress_phi,c='c', marker='o', linestyle='-.')
    plt.title('$\sigma_{\phi\phi}$ distribution')
    plt.xlabel('radius[mm]')
    plt.ylabel('$\sigma_{\phi\phi}$[MPa]')
    plt.show()
    
    
    plt.figure(4)
    plt.plot(time_fol,u_wide_1,c='b',label='node 1')
    plt.plot(time_fol,u_wide_4,c='m',label='node 4')
    plt.plot(time_fol,u_wide_7,c='r',label='node 7')
    plt.plot(time_fol,u_wide_11,c='g',label='$u_r(r=b$')
    plt.title('Time history of widening of the pipe')
    plt.xlabel('Time[s]')
    plt.ylabel('$u_r(r=b)$[mm]')
    plt.legend()
    plt.show()
    #plt.savefig('widening.png')
    
    plt.figure(5)
    plt.plot(r_elements,nodal_displacements,c='C0')
    plt.xlabel('radius[mm]')
    plt.ylabel('$u_r$[mm]')
    plt.title('Distribution of $u_r$')
    plt.axis([50,100,0.0425,0.066])
    plt.show()
    #plt.savefig('u_r.png')
    


