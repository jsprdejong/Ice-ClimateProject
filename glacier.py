from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np #This file is a preliminary version of the ice sheet model
import matplotlib.pyplot as plt

#set variables
Bs = 0 # Total surface mass budget
F = 0 # Calving flux (<0)
b0 = 3900 #m
s=0.1 # Bed slope
Lmax = b0/s
dx = 100. #m
beta=0.007 #m ice/a/m
alpha=3 #m^0.5
nu=10
E0=2900 # m height equilibrium line t=0
dt=1 #yr

ds_dl = 0

#set bottom slope glacier
x_arr = np.arange(0,Lmax,dx)

def mean_ice_thickness(alpha,nu,mean_s,L):
    return alpha/(1.+ nu*mean_s) * np.sqrt(L)

def bed_slope(x):
    return b0 - s * x

def dL_dt(L,alpha,nu,s,Bs,F):
    return (3*alpha/(2*(1+nu*s))*L**(1/2) \
    -alpha*nu/(1+nu*s)**2*L**(3/2)*ds_dl)**(-1)\
    *(Bs+F)        

# =============================================================================
# Main programm
# =============================================================================
b_arr = bed_slope(x_arr)
