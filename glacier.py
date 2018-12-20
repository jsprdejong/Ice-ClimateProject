from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np #This file is a preliminary version of the ice sheet model
import matplotlib.pyplot as plt
try: 
	from types import NoneType
except ImportError:
	print("Cannot import NoneType")
# =============================================================================
# Plot functions
# ============================================================================= 
def plot_results():
    plt.close("all")
    fig, ax = plt.subplots(nrows=1,ncols=2)
    fig.set_size_inches(12,6)
    ax[0].set_title("bed profile: "+str(bed_profile))
    ax[0].set_xlabel("x[km]")
    ax[0].set_ylabel("elevation[km]")
    
    b_arr = b(x_arr)
    ax[0].plot(x_arr/1000.,b_arr/1000.,"k")
    
    ax[0].set_xlim(0,max(x_arr/1000.))
    ax[0].set_ylim(0,max(b_arr/1000.)*1.1)

    ax[1].set_title("glacier length evolution, profile: "+str(bed_profile))
    ax[1].set_xlabel("time[yr]")
    ax[1].set_ylabel("length[km]")
    ax[1].plot(t_arr,L_arr/1000.)
    ax[1].set_xlim(0,max(t_arr))
    ax[1].set_ylim(0,max(L_arr/1000)*1.1)
    plt.show()

# =============================================================================
# Glacier model functions
# =============================================================================

def Hm(L):
    """mean glacier thickness"""
    mean_s_value = mean_s(L)
    return alpha/(1.+ nu*mean_s_value) * np.sqrt(L)

def b(x):
    """bed elevation"""
    if bed_profile == 'linear':
        return b0 - s * x
    elif bed_profile == 'concave':
        return ba + b0 * np.exp(-x/x_l)
    elif bed_profile == 'Aletschglacier':
        bed = np.zeros(len(x))
        bed[x < 8033] = b0 - s_A1 * x[x < 8033]
        bed[x >= 8033] = b1 - s_A2 * (x[x >= 8033] - 8033)
        return bed
  
def mean_b(L):
    if bed_profile == 'linear':
        return b0 - s*L/2.
    elif bed_profile == 'concave':
        return ba + x_l*b0/L*(1-np.exp(-L/x_l))
    elif bed_profile == 'Aletschglacier':
        if L < 8033:
            return b0 - s_A1 * L / 2.
        else:
            return (8033 / L) * b0 - s_A1 * 8033 / 2. + (1 - 8033 / L) * b1 - s_A2 * (L - 8033) / 2.
    	
def mean_s(L):
    """mean bed slope"""
    if bed_profile == 'linear':
        return s 
    if bed_profile == 'concave':
        return b0*(1-np.exp(-L/x_l))/L
    elif bed_profile == 'Aletschglacier':
        if L<8033:
            return s_A1
        else:
            return (8033.*s_A1+(L-8033.)*s_A2)/L

def dmean_s_dL(L,x):
    """change of the mean bed slope with respect to glacier length"""
    if bed_profile == 'linear':
        return 0
    elif bed_profile == 'concave':
        return -b0*(1-np.exp(-L/x_l))/L**2 + b0*x_l**(-1)/L
    elif bed_profile == 'Aletschglacier':
        if L < 8033:
            return 0
        else:
            return 8033.*(s_A2-s_A1)/L**2

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
    return (3. * alpha / (2. * (1. + nu * mean_s_value)) * L**(1./2.) \
    - alpha * nu / (1. + nu * mean_s_value)**2. * L**(3./2.) * ds_dl)**(-1.)\
    * (Bs + F)
  
# =============================================================================
# Glacier Tools
# =============================================================================
def integrate(L_prev,Hm_prev,Bs_prev,F_prev,E):
    """Integrate one timestep"""
    L_new = L_prev + dL_dt(L_prev,Bs_prev,F_prev)*dt
    Hm_new = Hm(L_new)
    Bs_new = Bs(Hm_new,E,L_new)
    F_new = F(Hm_new)    
    return L_new, Hm_new, Bs_new, F_new
  
def steady_state(E, L_0=0.001, min_dL_dt=0.1):
    """Returns year when steady-state is reached.
    
    :param E:          equilibrium line height
    :param L_0:        initial glacier length (default=0.001)
    :param min_dL_dt:  goal in difference in glacier length to 
                       achieve equilibrium
    """
    Bs_0 = 0
    Hm_0 = 0
    F_0= 0
    i = 0
    L_new, Hm_new, Bs_new, F_new = integrate(L_0,Hm_0,Bs_0,F_0,E)
    while dL_dt(L_new,Bs_new,F_new) > min_dL_dt or i<10:
        L_prev = L_new
        Hm_prev = Hm_new
        Bs_prev = Bs_new
        F_prev = F_new
        L_new, Hm_new, Bs_new, F_new = integrate(L_prev,Hm_prev,Bs_prev,F_prev,E)
        i+=1
    return i, L_new
  
def efolding(E,L_0):
    """Returns e-folding timescale.
    
    :param E: equilibrium line height
    :param L_0: initial glacier length
    """
    t_ss, L_ss = steady_state(E,L_0)
    L_efold = (1-1/np.exp(1))*(L_ss-L_0)+L_0

    Bs_arr = np.zeros(t_ss)
    Hm_arr = np.zeros(t_ss)
    F_arr = np.zeros(t_ss)
    L_arr = np.zeros(t_ss)
    L_arr[0] = L_0
    for i in range(t_ss-1):
        L_arr[i+1], Hm_arr[i+1], Bs_arr[i+1], F_arr[i+1] = \
        integrate(L_arr[i],Hm_arr[i],Bs_arr[i],F_arr[i])  
    L_diff = abs(L_arr-L_efold)

    t_efold = np.where(L_diff==min(L_diff))[0][0]*dt
    return t_efold
  
def find_current_ELA(L_today, dE = 50.):
    """This functions computes the current equlibirum height (ELA) 
    based on the current glacier length
    
    :param L_today: current glacier length
    :param dE: the pertubation of the ELA used in the newtons method (finite difference)
    """
    # initialize the ELA that is going to be found
    E_current = E0
    
    L_diff1 = np.inf
    while L_diff1 >10.:
        # steady state for current ELA
        year_steady1, L_steady1 = steady_state(E=E_current)
        L_diff1 = L_steady1-L_today
        
        # steady state for perturbated ELA
        year_steady2, L_steady2 = steady_state(E=E_current+dE)
        L_diff2 = L_steady2-L_today
        
        # derivative of the change
        deriv_L_diff = (L_diff2 - L_diff1)/dE
        
        # newton method 
        correction =  L_diff1/deriv_L_diff
        E_current = E_current - correction
  	
    return E_current
  
def project_ELA(year,dT_dt,current_ELA=None):
    """Projecting the evolution of the ELA
    
    :param year: year 
    :param dT_dt: the change of temperature per year
    """   
    if type(current_ELA) == NoneType:
        raise ValueError("The current ELA needs to be set or  computed!")
    return current_ELA + dE_dT *  dT_dt * (year - base_year)
  
def read_ELA():
    data = np.load("ELA.txt")
    years = data[:,0]
    ELA_perturbation = data[:,1]
    return years, ELA_perturbation

# =============================================================================
# Constants 
# =============================================================================
rho_w = 997. #kg/m3 density water
rho_i = 917. #kg/m3 density ice 

# =============================================================================
# Settings for the glacier model
# =============================================================================

# Aletschgletscher data
base_year = 2014 # year of measurement of the glacier length
L_2014 = 23000 # Length in 2014 (m)
E0 = 2900. # height equilibrium line t=0 (m)
W = 1. # m glacier width

# bed profile
b0 = 3900. # upper bound bed elevation (m) for linear case, upper bound for concave case would be b0+ba
ba = 0. # lower bound for concave bed profile
s = 0.1 # Bed slope
s_A1 = 0.1476401 #Bed slope upper part Aletschglacier
s_A2 = 0.0878518 #Bed slope lower part Aletschgacier
b1 = b0 - s_A1 * 8033 #bed elevation at slope division Aletschglacier
nu = 10.
alpha = 3. #m^0.5
beta = 0.007 #m ice/a/m
Lmax = b0/s
x_l = 7000. # m length scale for concave bed profile

# options for calving glaciers
d = 0 # water depth
kappa = 1. # fraction of mean ice thickness 
c = 1. # 1/yr calving parameter 

# ===================
# settings for the run
# ===================
dx = 100 #m
dt = 1 #yr
tmax = 100 #yr
bed_profile = 'Aletschglacier' # use 'linear' or 'concave' or 'Aletschglacier'
ds_dl = 0.

#set bottom slope glacier
x_arr = np.arange(0,Lmax,dx)
t_arr = np.arange(0,tmax,dt)

# =================
# climate scenario
# =================
dT_dt = 0.01 # change of temparture per year K/a 
dE_dT = 110.  # change of the ELA per temperature change m/K

# =============================================================================
# Main programm
# =============================================================================
Bs_arr = np.zeros(tmax)
Hm_arr = np.zeros(tmax)
F_arr = np.zeros(tmax)
L_arr = np.zeros(tmax)

# constant equilibrium height
E = E0

# initialize a first glacier
L_arr[0] = 0.01 # L_2014

for i in range(tmax-1):
    L_arr[i+1], Hm_arr[i+1], Bs_arr[i+1], F_arr[i+1] = \
    integrate(L_arr[i],Hm_arr[i],Bs_arr[i],F_arr[i],E)  

plot_results()