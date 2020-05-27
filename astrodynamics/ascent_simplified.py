# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:51:19 2020

@author: lucat
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from BAT import get_a, get_Mach, get_g, get_p, get_T, get_rho, get_ceff, Lambda, Mprop, get_Fl, get_Fd, get_TWratio, massfractions
from Aero.aero_calcs import aerodynamics_coefficients

#=====================================================================================================================================================================================================================
# Mars properties=  Atmospheric composition: 95.32% CO2, 2.7% N2, 1.6% Ar, 0.13% O2, 0.08% CO  
#=====================================================================================================================================================================================================================
Req=3396.2*10**3                #[m] equatorial radius
R=3389.5*10**3                  #[m] volumetric mean radius
g0_mars=3.71                    #[m/s^2] surface gravity
mu=0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter  
J2=1960.45*10**(-6) 
sol=24.6597*60*60               #[s] length of a day
i_mars=25.19*np.pi/180          #[rad] inlination of Mars' equator
omega_mars=0.7088218*10**(-4)   #[rad/s] Martian angular velocity
M_atm=43.34                     #[g/mol] Mean molecular weight
Rgas=8314.4621/M_atm               #[J/(kg*K)] Martian air gas constant
gamma=1.37                      #[-] heat capacity of Martian air
#=====================================================================================================================================================================================================================
# Node's properties & phasing orbit & Mars base
#=====================================================================================================================================================================================================================
i_node=42.5*np.pi/180           #[rad]
h_node=500*10**3                #[m]
i_phasing=i_node                #[rad]
V_phasing=3272.466              #[m/s]
h_phasing=609.74*10**3          #[m]
h0=-3*10**3                     #[m] altitude of the launch site wrt volumetric mean altitude
i0=i_node                       #[rad]
p0=get_p(h0)                        #[Pa]
T0=get_T(h0)                        #[K]
rho0=get_rho(p0,T0,Rgas)            #[kg/m^3]
g0=get_g(mu,Req,R,h0,i0,J2)      #[m/s^2]
#=====================================================================================================================================================================================================================
# Vehicle's properties
#=====================================================================================================================================================================================================================
Mdry=32251                     #[kg]
Isp=400                        #[s] LOX-LCH4
ceff=get_ceff(Isp)                 #[m/s]
Cl=0
Cd=0
d=7.67                         #[m] diameter of the vehicle assuming cylinder
S=np.pi/4*d**2                 #[m^2] reference surface area used here is the area of the circle of the cross-section
Ft=1000*10**3                  #[N] initial thrust at lift off to obtain T/W ratio=1.5
mdot0=881.5                    #[kg/s] 
tb=187                         #[s]


#=====================================================================================================================================================================================================================
# Simplified 2D point-mass model: gravity turn
#=====================================================================================================================================================================================================================

def dVxdt(T,D,L,M,Vz,Vx):
    dVxdt=(T-D)/M*Vx/np.sqrt(Vx*Vx+Vz*Vz)-L/M*Vz/np.sqrt(Vx*Vx+Vz*Vz)
    return dVxdt

#initial condition Vx: Vx(t=0):
Vx0=omega_mars*(R+h0)*np.cos(i0)

#final condition for Vx: Vx(t=te):
Vxe=V_phasing

def dVzdt(T,D,L,M,Vz,Vx,g):
    dVzdt=(T-D)/M*Vz/np.sqrt(Vx*Vx+Vz*Vz)-g+L/M*Vx/np.sqrt(Vx*Vx+Vz*Vz)
    return dVzdt

#initial condition Vz: Vz(t=0):
Vz0=0

#final condition for Vz: Vz(t=te):
Vze=0

#time span considered
dt=0.1


#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
Mp, Mwet=massfractions(ceff,198948,3272.466+916.29,45.842,289.485,88,30,391.492,373.65,487.238)
M=np.array([Mwet[0]])

Vx=np.array([Vx0])
Vz=np.array([Vz0])

Z=np.array([h0])
X=np.array([0])

p=np.array([get_p(Z[0])])
T=np.array([get_T(Z[0])])


while Z[-1]<h_phasing:
    t=t+dt
    
    #solve for Vx and Vz
    #compute V
    
    M=np.append(M, M[-1] - mdot0*dt)
    
    

#solving ODE's (Forward Euler) and updating new values
for i in range(int(tf/dt)):
    pnew=get_p(Z[i])
    p=np.append(p,pnew)
    Tnew=get_T(Z[i])
    T=np.append(T,Tnew)
    rho=get_rho(p,T,Rgas)
    a=get_a(gamma,Rgas,T)
    g=get_g(mu,Req,R,Z,i0,J2)
    TWratio=get_TWratio(Ft,M,g)
    Fd=get_Fd(Cd,rho,S,V)
    Fl=get_Fl(Cl,rho,S,V)
    
    #solve for Vx
    Vxnew=Vx[i]+dt*dVxdt(T,Fd[i],Fl[i],M[i],Vz[i],Vx[i])
    #solve for Vy
    Vznew=Vz[i]+dt*dVzdt(T,Fd[i],Fl[i],M[i],Vz[i],Vx[i],g[i])
    
    #update velocities
    Vx=np.append(Vx[i],Vxnew)
    Vz=np.append(Vz[i],Vznew)
    
    #absolute velocity
    V=np.sqrt(Vx*Vx+Vz*Vz)
    #Mach number
    Mach=get_Mach(V,a)
    #integrate the V's to get Xnew and Znew (note: h=Z)
    #Integrate using Trapezoidal rule 
    Znew=(t[i+1]-t[i])/2*(Vz[i]+Vznew)
    Z=np.append(Z,Znew)
    
    Xnew=(t[i+1]-t[i])/2*(Vx[i]+Vxnew)
    X=np.append(X,Xnew)
    
    #differentiate the V's to get ax and ay
    
    #compute new mass:
    Mnew=M[i]-mdot0
    #update M: 
    M=np.append(M,Mnew)
