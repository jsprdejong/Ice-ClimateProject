from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np #This file is a preliminary version of the ice sheet model
import matplotlib.pyplot as plt
# =============================================================================
# Plot functions
# ============================================================================= 
def plot_results():
    plt.close("all")
    fig, ax = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(12, 6)
    ax[0].set_title("bed profile: " + str(bed_profile))
    ax[0].set_xlabel("x[km]")
    ax[0].set_ylabel("elevation[km]")
    
    b_arr = b(x_arr)
    ax[0].plot(x_arr/1000., b_arr/1000., "k")
    
    ax[0].set_xlim(0, max(x_arr/1000.))
    ax[0].set_ylim(0, max(b_arr/1000.) * 1.1)

    ax[1].set_title("glacier length evolution, profile: "+str(bed_profile))
    ax[1].set_xlabel("time[yr]")
    ax[1].set_ylabel("length[km]")
    ax[1].plot(t_arr, L_arr/1000.)
    ax[1].set_xlim(0, max(t_arr))
    ax[1].set_ylim(0, max(L_arr/1000) * 1.1)
    plt.show()

# =============================================================================
# Glacier model functions
# =============================================================================

def Hm(L):
    """mean glacier thickness"""
    mean_s_value = mean_s(L)
    return alpha/(1. + nu*mean_s_value) * np.sqrt(L)

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
    
def balance_rate(x,E,Hm):
    return beta*(b(x) - E + Hm)
  
def mean_b(L):
    if bed_profile == 'linear':
        return b0 - s * L/2.
    elif bed_profile == 'concave':
        return ba + x_l * b0/L*(1-np.exp(-L/x_l))
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
            return (8033.*s_A1 + (L-8033.)*s_A2) / L

def dmean_s_dL(L):
    """change of the mean bed slope with respect to glacier length"""
    if bed_profile == 'linear':
        return 0
    elif bed_profile == 'concave':
        return -b0*(1-np.exp(-L/x_l))/L**2 + b0*x_l**(-1)/L*np.exp(-L/x_l)
    
    elif bed_profile == 'Aletschglacier':
        return np.heaviside(L-8033.,0)*8033.*(s_A2-s_A1)/L**2


def Bs(Hm, E, L):
    """Total surface mass budget"""  

    x = np.linspace(0,L,100)[1:]
    Delta_x = x[2] - x[1]
    
    W_arr = W_function(x)
    
    dot_b_arr = balance_rate(x,E,Hm)
    
    Delta_Bs_arr = np.zeros(len(W_arr))
    
    Delta_Bs_arr[:] = W_arr[:] * dot_b_arr[:] * Delta_x 
    
    Bs_sum = np.sum(Delta_Bs_arr)
    
    # analytical function
    #return beta * ( mean_b(L) + Hm - E) * 1 * L
    #numerical integration
    
    # Buckets (Jungfraufirn, Ewigschneefäld, Grüneggfirn respectively)
    return Bs_sum \
           + B_tribute(E, 6000., 1200, 3500, 2700, 3600) \
           + B_tribute(E, 4800., 1400, 3500, 2700, 3600) \
           + B_tribute(E, 2400., 350,  1500, 2700, 3600)

def B_tribute(E, L_trib, w_0_trib, w_end_trib, h_0_trib, h_end_trib):
    s_trib = (h_end_trib - h_0_trib)/L_trib
    q_trib = (w_end_trib - w_0_trib)/L_trib
    return beta * (w_0_trib * (h_0_trib - E) * L_trib + 1/2 * (s_trib * w_0_trib 
                   + (h_0_trib - E) * q_trib) * L_trib**2 + 1/3*s_trib*q_trib*L_trib**3)
    
def F(Hm,L):
    """"calving"""
    Hf= max(kappa*Hm, d*rho_w/rho_i) #ice thickness at glacier front
    W = Wm(L)
    return -c*d*Hf*W
    
def dL_dt(L,Bs,F):
    """prognostic equation for glacier length"""
    return (Bs + F)/Psi(L)
    
def Psi(L):
    mean_s_value = mean_s(L)
    dmean_s_dL_value = dmean_s_dL(L)
    
    Wm_value = Wm(L)
    dWm_dt_coef_value = dWm_dt_coef(L)
    A = 3. * alpha /(2. * (1. + nu * mean_s_value)) * L**(1./2.) * Wm_value
    B = alpha/(1+nu*mean_s_value)*L**(3/2)* dWm_dt_coef_value
    C = - alpha * nu / (1. + nu * mean_s_value)**2. * L**(3./2.) * Wm_value * dmean_s_dL_value
    
    return A+B+C

def W_function(x):
    return w0 + w1 * x *np.exp(-a*x)

def Wm(L):
    """ glacier width function"""
    return w0 + w1*a**(-2)*(L**(-1) - a* np.exp(-a*L) -L**(-1)*np.exp(-a*L))

def dWm_dt_coef(L):
    """returns the coefficent C in dWm/dt = C dL/dt"""
    return w1 * (- a**(-2)*L**(-2) + np.exp(-a*L) + a**(-2) * L**(-2) *np.exp(-a*L) + a**(-1)*L**(-1)*np.exp(-a*L))

# =============================================================================
# Glacier Tools
# =============================================================================
def integrate(L_prev, Hm_prev, Bs_prev, F_prev,E):
    """Integrate one timestep"""
    L_new = L_prev + dL_dt(L_prev,Bs_prev,F_prev)*dt
    Hm_new = Hm(L_new)
    Bs_new = Bs(Hm_new,E,L_new)
    F_new = F(Hm_new,L_new)    
    return L_new, Hm_new, Bs_new, F_new

def calc_glacier(E_guess):
    """ Compute evolution of the glacier.
    Returns Bs_arr, Hm_arr, F_arr, L_arr
    
    :param tmax:       integration time
    :E:                initial equilibrium height (scalar or array)
    :dT_dt:            temperature change per year
    :dE_dT:            change in equilibrium height per temperature change
    :L_0:              initial glacier length
    """
    
    # read ELA perturbations
    years, dELA = read_ELA()
    n_years = len(years)
    
    E_arr = np.zeros(n_years)
    E_arr[:] = E_guess + dELA[:]
    
    Bs_arr = np.zeros(n_years)
    Hm_arr = np.zeros(n_years)
    F_arr = np.zeros(n_years)
    L_arr = np.zeros(n_years)
    
    L_arr[0] = steady_state(E_guess) 
    
    for i in range(n_years-1):
        L_arr[i+1], Hm_arr[i+1], Bs_arr[i+1], F_arr[i+1] = \
        integrate(L_arr[i],Hm_arr[i],Bs_arr[i],F_arr[i],E_arr[i])
    return L_arr, years

def find_real_E0(E_guess, L_target, dE = 0.1, nu = 0.005):
    L_error1 = 10e5
    E_new = E_guess
    
    while abs(L_error1)>10.:
        
        E_guess = E_new
        
        L1arr, dummy = calc_glacier(E_guess)
        L2arr, dummy = calc_glacier(E_guess+dE)
        
        L1 = L1arr[-1]
        L2 = L2arr[-1]
        
        L_error1 = abs(L1 - L_target)
        L_error2 = abs(L2 - L_target)
        
        derror1_dE = (L_error2-L_error1)/dE
        
        E_new = E_guess - nu * L1 / derror1_dE
    
    return E_guess

def project_future_glacier(tmax, E, dT_dt=0, dE_dT=0, L_0=10.01, \
                           year_start = False, year_stab = False):
    """ Compute evolution of the glacier.
    Returns Bs_arr, Hm_arr, F_arr, L_arr
    
    :param tmax:       integration time
    :E:                initial equilibrium height
    :dT_dt:            temperature change per year
    :dE_dT:            change in equilibrium height per temperature change
    :L_0:              initial glacier length
    """
    
    if year_start != False:
        years = np.arange(year_start, year_start + tmax,1)
    else:
        years = np.arange(0, tmax, 1)
    
    Bs_arr = np.zeros(tmax)
    Hm_arr = np.zeros(tmax)
    F_arr = np.zeros(tmax)
    L_arr = np.zeros(tmax)

    L_arr[0] = L_0 # L_2014

    for i in range(tmax - 1):
        if year_stab!=False:
            if years[i] == year_stab:
                dT_dt = 0
        L_arr[i+1], Hm_arr[i+1], Bs_arr[i+1], F_arr[i+1] = \
        integrate(L_arr[i],Hm_arr[i],Bs_arr[i],F_arr[i],E)
        E = E + dE_dT * dT_dt
    
    return L_arr, years

def steady_state(E, L_0=10.001, Delta_L = 10., nu=0.1):
    
    Hm_0 = Hm(L_0)
    Bs_0 = Bs(Hm_0,E,L_0)
    F_0 = F(Hm_0,L_0) 
    
    for i in range(300):
        L_0, Hm_0, Bs_0, F_0 = integrate(L_0, Hm_0, Bs_0, F_0, E)
    
    correction = 1000.
    while abs(correction) > 10.:
        Hm_0 = Hm(L_0)

        Bs_0 = Bs(Hm_0,E,L_0)
        F_0 = F(Hm_0,L_0) 
        
        L_1 = L_0 + Delta_L
        Hm_1 = Hm(L_1)
        Bs_1 = Bs(Hm_1,E,L_1)
        F_1 = F(Hm_1,L_1)  
    
        dL_dt0 = dL_dt(L_0, Bs_0, F_0)
        dL_dt1 = dL_dt(L_1, Bs_1, F_1)
        
        ddL_dtdL = (dL_dt1 - dL_dt0)/Delta_L 
        
        correction = dL_dt0/ddL_dtdL
        L_0 = L_0 - nu * correction
        if L_0 <= 0.:
            L_0 = 0
            break
    return L_0
  
def efolding(E, L_ref):
    """Returns e-folding timescale for the adjustment to a new equilibium height (E).
    
    :param E: new equilibrium line height
    :param L_ref: glacier length for an reference equilibrium height
    """
    L_ss = steady_state(E, L_ref)
    
    if abs(L_ss - L_ref) < 10:
        return np.nan
    
    L_efold =  (1 - np.exp(-1)) * (L_ss - L_ref) + L_ref

    L_new,L_old = L_ref,L_ref
    H_new, B_new, F_new = 0, 0, 0
    
    t_efold = 0
    
    while (L_new > L_efold and L_old > L_efold) or (L_new < L_efold and L_old < L_efold):
        L_old, H_old, B_old, F_old = L_new, H_new, B_new, F_new 
        L_new, H_new, B_new, F_new = integrate(L_old, H_old, B_old, F_old,E)  
        
        t_efold+=1
        
    return t_efold
  
def read_ELA():
    data = np.loadtxt("ELA.txt")
    years = data[:-4,0]
    ELA_perturbation = data[:-4,1]
    return years, ELA_perturbation

def E_fixed_points(L_0=10.001):
  """ Returns equilibrium glacier length and integration time for
      different values of E
  """
  E_arr = np.arange(1000, 4000, 10)
  lmax_arr = np.array([steady_state(E,L_0) for E in E_arr])
  return E_arr, lmax_arr

def E_vs_efolding(ref_E):
  """ Calculate efolding timescales for different values of E when 
      the glacier length is the equilibrium length calculated
      from a reference value of E (E0)
  """
  L_equil = steady_state(ref_E)
  E_arr = np.arange(1000, 4000, 10)
  t_efold_arr = np.array([efolding(E,L_equil) for E in E_arr])
  return E_arr, t_efold_arr

#%% =============================================================================
# Constants 
# =============================================================================
rho_w = 997. #kg/m3 density water
rho_i = 917. #kg/m3 density ice 

# =============================================================================
# Settings for the glacier model
# =============================================================================
# Aletschgletscher data
base_year = 2014 # year of measurement of the glacier length
L_2014 = 23950 # Length in 2014 (m)
E0 = 2900. # height equilibrium line t=0 (m)
w0 = 1800.
w1 = 0 # 3.
a = 0.00045


# bed profile
b0 = 3900. # upper bound bed elevation (m) for linear case, upper bound for concave case would be b0+ba
ba = 0. # lower bound for concave bed profile
s = 0.1#476401  # Bed slope
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
plt.figure(7)
plt.plot(x_arr,W_function(x_arr))

# initialize a first glacier
L_arr, years = project_future_glacier(tmax,E0,L_0=10.1)
plot_results()

#%%
plt.figure(2)
plt.title("E-folding time ")
plt.xlabel("E [m]")
plt.ylabel("e-folding time [yr]")
E_arr, t_efold_arr = E_vs_efolding(E0)
plt.plot(E_arr, t_efold_arr)

plt.savefig("figures/efolding.png",dpi=300)




plt.figure(3)
plt.title("Stable steady state")
plt.xlabel("E [m]")
plt.ylabel("L [km]")
E_arr, lmax_arr = E_fixed_points()
plt.plot(E_arr, lmax_arr/1000)

plt.savefig("figures/steady_states.png",dpi=300)





E0 = find_real_E0(2900.,21500.)
L_arr, year_arr = calc_glacier(E0)

plt.figure(4)
plt.plot(year_arr,L_arr)
plt.title("Historical Glacier Length Evolution")
plt.xlabel("t [yr]")
plt.ylabel("L [km]")
plt.savefig("figures/historical_glacier.png",dpi=300)






E_2014 = E0 + read_ELA()[1][-1]
y_start = 2014
y_end = 2200
y_stab = 2100
L_arr_future_001, years_future = \
project_future_glacier(y_end-y_start, E_2014, 0.01, 110, L_2014, y_start, y_stab)
L_arr_future_002, years_future = \
project_future_glacier(y_end-y_start, E_2014, 0.02, 110, L_2014, y_start, y_stab)
L_arr_future_004, years_future = \
project_future_glacier(y_end-y_start, E_2014, 0.04, 110, L_2014, y_start, y_stab)

plt.figure(5)
plt.title("Glacier Length Evolution under Climate Change Scenarios")
plt.xlabel("t [yr]")
plt.ylabel("L [km]")
plt.axvline(y_stab, linestyle='--', c='red', label='stabilization year')
plt.plot(years_future, L_arr_future_001, label='0.01 K/a')
plt.plot(years_future, L_arr_future_002, label='0.02 K/a')
plt.plot(years_future, L_arr_future_004, label='0.04 K/a')
plt.legend()
plt.savefig("figures/future_glacier.png", dpi=300)






ELA_years, ELA_perturbation = read_ELA()
plt.figure(6)
plt.title("Equilibrium Line Altitude Data")
plt.xlabel("t [yr]")
plt.ylabel("ELA perturbation (m)")
plt.plot(ELA_years, ELA_perturbation)
plt.savefig("figures/ELA_evolution.png",dpi=300)
