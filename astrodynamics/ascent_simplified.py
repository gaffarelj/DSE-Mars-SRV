# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:51:19 2020

@author: lucat
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from BAT import get_a, get_Mach, get_g, get_p, get_T, get_rho, get_ceff, Lambda, Mprop, get_Fl, get_Fd, get_TWratio, massfractions
#from Aero.aero_calcs import aerodynamics_coefficients
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
from aero_calcs import aerodynamics_coefficients
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
Rgas=8314.4621/M_atm            #[J/(kg*K)] Martian air gas constant
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
    dVxdt=(T-D)/M*Vx/(np.sqrt(Vx*Vx+Vz*Vz))-L/M*Vz/(np.sqrt(Vx*Vx+Vz*Vz))
    return dVxdt

#free Vx due to the rotation of Mars: Vx(t=0):
Vxfree=omega_mars*(R+h0)*np.cos(i0)

#initial condition for Vx wrt to surface: Vx(t=0):
Vx0=0.00001

#final condition for Vx: Vx(t=te):
Vxe=V_phasing

def dVzdt(T,D,L,M,Vz,Vx,g):
    dVzdt=(T-D)/M*Vz/(np.sqrt(Vx*Vx+Vz*Vz))-g+L/M*Vx/(np.sqrt(Vx*Vx+Vz*Vz))
    return dVzdt

#initial condition Vz: Vz(t=0):
Vz0=0.00001

#final condition for Vz: Vz(t=te):
Vze=0

#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
Mp, Mwet=massfractions(ceff,198948,3272.466+916.29,45.842,289.485,88,30,391.492,373.65,487.238)
M=np.array([Mwet[0]])

Vx=np.array([Vx0])
Vz=np.array([Vz0])
V=np.sqrt(Vx**2+Vz**2)

#earth reference frame is on the surface of Mars (where the base is, i.e. -3km)
Z=np.array([0])
X=np.array([0])


dt=0.01
t_tot=0

Mach_list=[]
ax_list=[]
az_list=[]


while Z[-1]<h_phasing+abs(h0):
#for i in range(1):
    t_tot=t_tot+dt #k
    p=get_p(Z[-1]) #k
    T=get_T(Z[-1]) #k
    rho=get_rho(p,T,Rgas) #k
    a=get_a(gamma,Rgas, T) #k
    Mach=get_Mach(V[-1],a) #k
    Mach_list.append(Mach)
    Cl,Cd=aerodynamics_coefficients(Mach,0)
    g=get_g(mu,Req,R,Z[-1],i0,J2)
    Fd=get_Fd(Cd,rho,S,V[-1])
    Fl=get_Fl(Cl,rho,S,V[-1])
    
    
    #solve for Vx
    ax=dVxdt(Ft,Fd,Fl,M[-1],Vz[-1],Vx[-1])
    ax_list.append(ax)
    Vxnew=Vx[-1]+dt*ax
    #Integrate using Trapezoidal rule 
    Xnew=dt/2*(Vx[-1]+Vxnew)
    X=np.append(X,Xnew)
    
    #solve for Vz
    az=dVzdt(Ft,Fd,Fl,M[-1],Vz[-1],Vx[-1],g)
    az_list.append(az)
    Vznew=Vz[-1]+dt*az
    #Integrate using Trapezoidal rule 
    Znew=dt/2*(Vz[-1]+Vznew)
    Z=np.append(Z,Znew)
    print(az)

    #update velocities
    Vx=np.append(Vx,Vxnew)
    Vz=np.append(Vz,Vznew)
    #compute V
    V=np.append(V,np.sqrt(Vxnew**2+Vznew**2))
    
  
    M=np.append(M, M[-1] - mdot0*dt)

"""
t_array=np.linspace(0,t_tot,t_tot/dt)  
plt.plot(t_array,V)    
plt.show()
"""  
