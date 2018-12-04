from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np #This file is a preliminary version of the ice sheet model
import matplotlib.pyplot as plt

#set variables
b0 = 3900 #m
s=0.1
Lmax = b0/s
dx = 100. #m
beta=0.007 #m ice/a/m
alpha=3 #m^0.5
nu=10
E0=2900 # m height equilibrium line t=0
dt=1 #yr


#set bottom slope glacier
x_arr = np.arange(0,Lmax,dx)

def bed_slope(x):
    return b0 - s * x

# =============================================================================
# Main programm
# =============================================================================
b_arr = bed_slope(x_arr)

print(dt)