from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np #This file is a preliminary version of the ice sheet model
import matplotlib.pyplot as plt

#set constants
b0 = 3900. # m
ba = 0 #m
s = 0.1 # Bed slope
Lmax = b0/s
beta=0.007 #m ice/a/m
alpha=3. #m^0.5
nu=10.
E0=2900. # m height equilibrium line t=0
rho_w = 997. #kg/m3 density water
rho_i = 917. #kg/m3 density ice 
c=1. # 1/yr calving parameter 
x_l = 7000. #m length scale where bed profile drops off
W=1. #m glacier width
d = 0# water depth
kappa = 1. # fraction of mean ice thickness 


#settings
dx = 100 #m
dt=1 #yr
tmax = 100 #yr
bed_profile = 'linear' # use 'linear' or 'concave'
ds_dl = 0.

#set bottom slope glacier
x_arr = np.arange(0,Lmax,dx)

def Hm(L):
    """mean glacier thickness"""
    mean_s_value = mean_s(L)
    return alpha/(1.+ nu*mean_s_value) * np.sqrt(L)

def b(x):
    """bed elevation"""
    return b0 - s * x
  
def mean_b(L):
    if bed_profile == 'linear':
  		return b0 - s*L/2.
    elif bed_profile == 'concave':
        return ba + x_l*b0/L*(1-np.exp(-L/x_l))
    
def mean_s(L):
    """mean bed slope"""
    if bed_profile == 'linear':
  		return s 
    if bed_profile == 'concave':
        	return b0(1-np.exp(-L/x_l))/L

def dmean_s_dL(L,x):
    """change of the mean bed slope with respect to glacier length"""
    if bed_profile == 'linear':
        return 0
    elif bed_profile == 'concave':
        return -b0(1-np.exp(-L/x_l))/L**2 + b0*x_l**(-1)/L

def Bs(Hm,E,L):
    """Total surface mass budget"""  
    mean_b_value = mean_b(L)
    return beta * (mean_b_value + Hm - E) * W * L
  
def F(Hm):
    """"calving"""
    Hf= max(kappa*Hm, d*rho_w/rho_i) #ice thickness at glacier front
    return -c*d*Hf*W

def dL_dt(L,Bs,F):
    """prognostic equation for glacier length"""
    mean_s_value = mean_s(L)  
    return (3 * alpha / (2 * (1 + nu * mean_s_value)) * L**(1/2) \
    - alpha * nu / (1 + nu * mean_s_value)**2 * L**(3/2) * ds_dl)**(-1)\
    * (Bs + F)   

  
# =============================================================================
# Main programm
# =============================================================================

Bs_arr = np.zeros(tmax)
Hm_arr = np.zeros(tmax)
F_arr = np.zeros(tmax)
L_arr = np.zeros(tmax)

E = E0
for i in range(tmax-1):
    L_arr[i+1] = L_arr[i] + dL_dt(L_arr[i],Bs_arr[i],F_arr[i])*dt
    Hm_arr[i+1] = Hm(L_arr[i+1])
    Bs_arr[i+1] = Bs(Hm_arr[i+1],E,L_arr[i+1])
    F_arr[i+1] =  F(Hm_arr[i+1])    

