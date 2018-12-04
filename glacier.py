import numpy as np #This file is a preliminary version of the ice sheet model
import matplotlib.pyplot as plt

#set variables
b0=3900 #m
s=0.1
beta=0.007 #m ice/a/m
alpha=3 #m^0.5
nu=10
E0=2900 # m height equilibrium line t=0
dt=1 #yr

#set bottom slope glacier
x=np.zeros(len(b0/s))
b=np.zeros(len(x))

#bottom bed slope
for i in range(len(x)):
    b[i]=b0-sx