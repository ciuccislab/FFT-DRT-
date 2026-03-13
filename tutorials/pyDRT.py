#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 11:15:47 2022

@author: francesco
"""
import numpy as np
from math import pi, log, exp, sqrt
from scipy.optimize import fsolve
from scipy import integrate


def compute_A_re(freq_vec, tau_vec):
    
    omega_vec = 2.*pi*freq_vec
    N_f = freq_vec.size
    N_tau = tau_vec.size

    out_A_re = np.zeros((N_f, N_tau))
    Delta_log_tau_vec = np.zeros(N_tau)
    
    for n in range(0, N_tau):
        if n == 0:
            Delta_log_tau_vec[n] = log(tau_vec[n+1]/tau_vec[n])
        elif n == N_tau-1:
            Delta_log_tau_vec[n] = log(tau_vec[n]/tau_vec[n-1])
        else:
            Delta_log_tau_vec[n] = log(tau_vec[n+1]/tau_vec[n-1])

    for n in range(0, N_tau):      
        out_A_re[:,n] = 0.5/(1+(omega_vec*tau_vec[n])**2)*Delta_log_tau_vec[n]
            
    return out_A_re


def compute_A_im(freq_vec, tau_vec):
    
    omega_vec = 2.*pi*freq_vec

    N_tau = tau_vec.size
    N_f = freq_vec.size

    out_A_im = np.zeros((N_f, N_tau))
    Delta_log_tau_vec = np.zeros(N_tau)

    for n in range(0, N_tau):
        if n == 0:
            Delta_log_tau_vec[n] = log(tau_vec[n+1]/tau_vec[n])
        elif n == N_tau-1:
            Delta_log_tau_vec[n] = log(tau_vec[n]/tau_vec[n-1])
        else:
            Delta_log_tau_vec[n] = log(tau_vec[n+1]/tau_vec[n-1])

    for n in range(0, N_tau):            
        out_A_im[:,n] = -0.5*(omega_vec*tau_vec[n])/(1+(omega_vec*tau_vec[n])**2)*Delta_log_tau_vec[n]
            
    return out_A_im


def compute_L2(tau):
    
    N_tau = tau.size
    log_tau = np.log(tau)
    out_L2 = np.zeros((N_tau-2, N_tau))
    
    for n in range(0, N_tau-2):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L2[n,n] = 1./(delta_loc**2)
        out_L2[n,n+1] = -2./(delta_loc**2)
        out_L2[n,n+2] = 1./(delta_loc**2)


    return out_L2


def compute_L1(tau):
    
    N_tau = tau.size
    log_tau = np.log(tau)
    out_L1 = np.zeros((N_tau-1, N_tau))
    
    for n in range(0, N_tau-1):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L1[n,n] = 1./delta_loc
        out_L1[n,n+1] = -1./delta_loc
    
    return out_L1


def compute_L2_aristo(tau, sigma_aristo):
    
    N_tau = tau.size
    log_tau = np.log(tau)
    out_L2 = np.zeros((N_tau, N_tau))
    
    for n in range(0, N_tau-2):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L2[n+1,n] = 1./(delta_loc**2)
        out_L2[n+1,n+1] = -2./(delta_loc**2)
        out_L2[n+1,n+2] = 1./(delta_loc**2)
        
    # aristotelian
    out_L2[0, 0] = sigma_aristo
    out_L2[N_tau-1, N_tau-1] = sigma_aristo

    return out_L2


def compute_L1_aristo(tau, sigma_aristo):
    
    N_tau = tau.size
    log_tau = np.log(tau)
    out_L1 = np.zeros((N_tau+1, N_tau))
    
    for n in range(0, N_tau-1):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L1[n+1,n] = 1./delta_loc
        out_L1[n+1,n+1] = -1./delta_loc

    # aristotelian
    out_L1[0, 0] = sigma_aristo
    out_L1[N_tau, N_tau-1] = sigma_aristo
    
    return out_L1


def compute_P_q_norm(A, Z, freq_plus):
    
    A_tilde = np.zeros_like(A)
    N_f = freq_plus.shape[0]
    for m, freq_loc in enumerate(freq_plus):
        A_tilde[m,:] = 1.0/freq_loc*A[m,:]
        A_tilde[m+N_f,:] = 1.0/freq_loc*A[m+N_f,:]

    out_P = A.T@A_tilde
    out_q = -Z.T@A_tilde

    return out_P, out_q


def compute_P_q(A, Z, freq_plus):
        
    out_P = A.T@A
    out_q = -Z.T@A

    return out_P, out_q


def mittag_leffler_func(x, phi): # fourth-order global Padé approximation of the Mittag-Leffler function
    
    # reference: Highly accurate global Padé approximations of generalized Mittag–Leffler function and its inverse,
    # I.O. Sarumi, K.M. Furati, A.Q.M. Khaliq, Journal of Scientific Computing. 82 (2020) 46.
    
    # define the 7x7 matrix Gamma_alpha,beta that contains only 23 non-zero entries (see (31) in the reference)
    Gamma_mat = np.zeros((7,7))
    
    for n in range(3):
        Gamma_mat[n,n] = 1
        Gamma_mat[n+3,n+3] = 1/gamma(4*phi+1)
    
    for n in range(0, 4):
        Gamma_mat[n,n+3] = -1/gamma(phi+1)
        Gamma_mat[n+1,n+3] = 1/gamma(2*phi+1)
        Gamma_mat[n+2,n+3] = -1/gamma(3*phi+1)
    
    Gamma_mat[4,3] = Gamma_mat[5,4] = -1/gamma(5*phi+1)
    Gamma_mat[5,3] = 1/gamma(6*phi+1)
    Gamma_mat[6,2] = 1
    Gamma_mat[6,6] = -1
    
    # compute the parameters p1, p2, p3, q0, q1, q2, and q3 (see (31) in the reference)
    y_pq = np.zeros(7)
    y_pq[3] = -1
    y_pq[4] = 1/gamma(phi+1)
    y_pq[5] = -1/gamma(2*phi+1)
    y_pq[6] = -1/gamma(1-phi)
    pq_vec = np.linalg.solve(Gamma_mat, y_pq)
    
    # define the Pade approximate of E_alpha,beta(-x) as per (26) in the reference
    num = pq_vec[0] + pq_vec[1]*x + pq_vec[2]*x**2 + x**3
    denom = pq_vec[3] + pq_vec[4]*x + pq_vec[5]*x**2 + pq_vec[6]*x**3 + x**4
    
    R_out = num/denom    
    
    return R_out


def rho_function(t, tau1, tau2, tau3, t0, delta_time, current_type): # see (12) and (16)
    
    if current_type == 'pulse': # pulse current
    
        rho_out = 0.5*(1-np.exp(-t/tau1))*np.log(tau2/tau3) 
        
    else: # triangular current
        
        if t0-delta_time<=t<=t0:
        
            rho_out = 0.5*(1 + (t-t0-tau1)/delta_time + tau1/delta_time*np.exp(-(t-t0+delta_time)/tau1))*np.log(tau2/tau3)
        
        elif t0-delta_time<=t<=t0:
        
            rho_out = 0.5*(1 + (t0-t+tau1)/delta_time - (1+tau1/delta_time)*np.exp(-(t-t0)/tau1))*np.log(tau2/tau3)
        
        else:
        
            rho_out = 0
    
    return rho_out


def assemble_B(time_vec, tau_vec, t0, delta_time, current_type): 
    
    # B_matrix defined in (10), (11), and (16) in the DRT-DRT theory
    
    M_pulse = len(time_vec)
    N_tau = len(tau_vec)
        
    out_B = np.zeros([M_pulse, N_tau+1])
    
    for m in range(M_pulse):
        
        out_B[m, 0] = 1 # see (10)
        
        for n in range(1, N_tau+1):
            
            if current_type == 'pulse': # pulse current (see (11))
                
                if time_vec[m] < delta_time:
                
                    if n==1:
                    
                        out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n], tau_vec[n-1], t0, delta_time, current_type)
                    
                    elif n==N_tau:
                    
                        out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n-1], tau_vec[n-2], t0, delta_time, current_type)

                    else:
                    
                        out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n], tau_vec[n-2], t0, delta_time, current_type)
            
                else:
                
                    if n==1:
                    
                        out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n], tau_vec[n-1], t0, delta_time, current_type) -\
                        rho_function(time_vec[m]-delta_time, tau_vec[n-1], tau_vec[n], tau_vec[n-1], t0, delta_time, current_type) 
                    
                    elif n==N_tau:
                    
                        out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n-1], tau_vec[n-2], t0, delta_time, current_type) -\
                        rho_function(time_vec[m]-delta_time, tau_vec[n-1], tau_vec[n-1], tau_vec[n-2], t0, delta_time, current_type) 
                    
                    else:
                    
                        out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n], tau_vec[n-2], t0, delta_time, current_type) -\
                        rho_function(time_vec[m]-delta_time, tau_vec[n-1], tau_vec[n], tau_vec[n-2], t0, delta_time, current_type) 
               
            else: # triangular current (see (16))
            
                if n==1:
                    
                    out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n], tau_vec[n-1], t0, delta_time, current_type)
                    
                elif n==N_tau:
                    
                    out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n-1], tau_vec[n-2], t0, delta_time, current_type)

                else:
                    
                    out_B[m, n] = rho_function(time_vec[m], tau_vec[n-1], tau_vec[n], tau_vec[n-2], t0, delta_time, current_type)
            
    return out_B