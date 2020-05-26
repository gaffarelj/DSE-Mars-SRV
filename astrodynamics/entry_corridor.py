# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:22:41 2020

@author: lucat
"""

import numpy as np
from matplotlib import pyplot as plt 
from mars_standard_atmosphere import *
from reentry import height, velocity


#input parameters
steps=10000
g0=3.71                 #[m/s]
h_final=80*10**3        #[m]
h_LMO=500*10**3         #[m]
mu=0.042828*1000000*(10**3)**3 #[m^3/s^2]
R_mars=3389.5*10**3         #[m] mean volumetric radius of Mars
a_LMO=R_mars+h_LMO
V_LMO=np.sqrt(mu/a_LMO)
R = 8314.4621/44.1
gamma = 1.37
p0=get_pressure(0)
T0=get_temperature(0)
rho0=get_density(p0, T0, R=R)
Cl=1.5 
Cd=2.75 
m=20000                 #[kg]
S=65      
n=4


h=np.linspace(0,h_final,steps)
V=np.linspace(0,3500,steps)

#considering only equilibrium-glide condition and g-load condition
#hence no stagnation heat flux condition
def rho_eq_glide(Veq,W,S,Cl,Vc):
    return 2*(W/S)/Cl*(1/(Veq*Veq)-1/(Vc*Vc))

def rho_g(Vg,n,m,g0,S,Cd,Cl):
    return 2*n*(m*9.81)/(Vg*Vg*S*np.sqrt(Cd*Cd+Cl*Cl))

def g(g0,h,R):
    return g0*1/(1+h/R)**2

def gamma_bar(Hs,R,Cl,Cd,Vc,V):
    return np.arcsin(-1/(1/Hs*R)*2/(Cl/Cd)*(Vc*Vc)/(V*V))

#exponential model for Mars atmosphere
def h_exp(Hs,rho,rho0):
    return -Hs*np.log(rho/rho0)
"""
Hs=[] 
for i in range(len(h)):
    Hs.append((R*get_temperature(h[i]))/g0)
Hs=np.array(Hs)
"""
Hs=8.8*10**3 #[m] 8.8 km->from Reentry reader, 11.1->from wikipedia

#Vc=np.sqrt(g(g0,h,R)*R) #not a valid assumption!
Vc=np.sqrt(mu/(R_mars+h))
W=m*g(g0,h,R)
rho_g=rho_g(V,n,m,g0,S,Cd,Cl)
rho_eq_glide=rho_eq_glide(V,W,S,Cl,Vc)


h_eq_glide=h_exp(Hs,rho_eq_glide,rho0)
h_g=h_exp(Hs,rho_g,rho0)

gamma_bar=gamma_bar(Hs,R,Cl,Cd,Vc,V)

#trajectories of different periapsis. Naming h_periapsis, where periapsis=altitude of periapsis computed for

#periapsis=400km
h_400=[31401,70000,10004]
V_400=[1113.88,3462.96,369.53]

#periapsis=300km
h_300=[36100,7000,10000]
V_300=[1430.04,3480.15,369.85]

plt.rcParams.update({'font.size': 12})
#plot
plt.plot(V,h_eq_glide,color="navy",label="skipping flight constraint")
plt.plot(V,h_g,color="firebrick",label="g-load constraint")
plt.fill_between(V, h_eq_glide, h_g,color="ivory")
#plt.plot(V_400,h_400,"r+")
plt.plot(velocity, height,"--",color="mediumseagreen",label="re-entry flight envelope: periapsis=3089 km")
plt.grid(color="gainsboro")
plt.ylim(bottom=0)
plt.title("Entry Corridor: Cl="+ str(Cl) +", Cd=" + str(Cd) + ", S=" + str(round(S,2)) + " [m^2], m=" + str(round(m,2)) + " [kg]")
plt.xlabel("velocity [m/s]")
plt.ylabel("altitude [m]")
plt.legend()
plt.show()

