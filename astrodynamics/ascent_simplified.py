# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:51:19 2020

@author: lucat
"""

import numpy as np
import matplotlib.pyplot as plt
from BAT import get_a, get_Mach, get_g, get_p, get_T, get_rho, get_ceff, Lambda, Mprop, get_Fl, get_Fd, get_Ft, get_TWratio, massfractions,ROM_deltaV
#from Aero.aero_calcs import aerodynamics_coefficients
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
from aero_calcs import H_aerodynamics_coefficients
import time
start_time = time.time()

#=====================================================================================================================================================================================================================
# Mars properties =  Atmospheric composition: 95.32% CO2, 2.7% N2, 1.6% Ar, 0.13% O2, 0.08% CO  
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
gamma_gas=1.37                      #[-] heat capacity of Martian air
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
Isp= 350                       #[s] LOX-LCH4
ceff=get_ceff(Isp)             #[m/s]
Cl=0
Cd=0
d=7.67                         #[m] diameter of the vehicle assuming cylinder
n=6                            #[-] number of engines
S=np.pi/4*d**2                 #[m^2] reference surface area used here is the area of the circle of the cross-section
De=2                           #[m] Diameter of an engine
Ae=np.pi/4*De*De               #[m] Exit area of engine
pe=5066.25                     #[Pa] exit pressure of exhaust gasses
                      


#=====================================================================================================================================================================================================================
# Switches
#=====================================================================================================================================================================================================================

plotting=True       #Do you wanna plot? no=False

#=====================================================================================================================================================================================================================
# Simplified 2D point-mass model: gravity turn
#=====================================================================================================================================================================================================================

def dVxdt(T,D,L,M,Vz,Vx):
    dVxdt=(T-D)/M*Vx/(np.sqrt(Vx*Vx+Vz*Vz))-L/M*Vz/(np.sqrt(Vx*Vx+Vz*Vz))
    return dVxdt

#free Vx due to the rotation of Mars: Vx(t=0):
Vxfree=omega_mars*(R+h0)*np.cos(i0)

#final condition for Vx: Vx(t=te):
Vxe=V_phasing

def dVzdt(T,D,L,M,Vz,Vx,g):
    dVzdt=(T-D)/M*Vz/(np.sqrt(Vx*Vx+Vz*Vz))-g+L/M*Vx/(np.sqrt(Vx*Vx+Vz*Vz))
    return dVzdt



#initial conditions
V0=0.00001                      #8.33[m/s] circa 30 km/h
gamma0=(90-3.2)*np.pi/180      #degs off vertical axis 3.2

Vx0=V0*np.cos(gamma0)
Vz0=V0*np.sin(gamma0)
#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
Mp, Mwet=massfractions(ceff,198948,3272.466+916.29,45.842,289.485,88,30,391.492,373.65,487.238)
M=np.array([Mwet[0]])
    
Vx=np.array([Vx0])
Vz=np.array([Vz0])
V=np.array([V0])
    
#earth reference frame is on the surface of Mars (where the base is, i.e. -3km)
Z=np.array([0.0])    #Z0=10m
X=np.array([0.0]) #X0=0.524m

TW0=1.5
TWe=4
    
g=np.array([get_g(mu,Req,R,Z[-1],i0,J2)])
dt=0.01
t_tot=0
    
 
p=np.array([get_p(Z[-1])])
T=np.array([get_T(Z[-1])])
rho=np.array([get_rho(p[-1],T[-1],Rgas)])
Mach=np.array([get_Mach(V[0],get_a(gamma_gas,Rgas, T[0]))])
#TWratio=np.array([get_TWprofile(mode,0,g[0],1.5,4,tb,i0,Z[0])])
tb=ceff/(TW0*9.80665/g0_mars*9.80665)*np.log(M[0]/(M[0]-Mp[0]))+41
TWratio=np.linspace(TW0*9.80665/g0_mars,TWe*9.80665/g0_mars,round(tb/dt))
Cl0,Cd0=H_aerodynamics_coefficients(Mach,0)
Cl=np.array([Cl0])
Cd=np.array([Cd0])
Fd=np.array([get_Fd(Cd[0],rho[0],S,V[0])])
Fl=np.array([get_Fl(Cl[0],rho[0],S,V[0])])
mdot=np.array([1/ceff*(M[0]*g[0]*TWratio[0]-Ae*(pe-p[0]))])
Ft=np.array([get_Ft(mdot[0],ceff,Ae,pe,p)])
ax_array=np.array([dVxdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1])])
az_array=np.array([dVzdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1],g[-1])])

#TWratio=np.array([Ft[0]/(M[0]*g[0])])
    
i=-1
    
while Z[-1]<h_phasing+abs(h0) and Z[-1]>=0:
#for i in range(1):
        
        #######################################################
        #   Compute parameters needed                         #
        #######################################################
        t_tot+=dt 
        i+=1
        
        a=get_a(gamma_gas,Rgas, T[-1]) 
        #######################################################
        #   Solve for ax, az, Vx, Vz, X, Z                    #
        #######################################################
        #solve for Vx
        ax=dVxdt(Ft[-1],Fd[-1],Fl[-1],M[-1],Vz[-1],Vx[-1])
        ax_array=np.append(ax_array,ax)
        Vxnew=Vx[-1]+dt*ax
        
        #Integrate using Trapezoidal rule 
        deltaX=dt/2*(Vx[-1]+Vxnew)
        Xnew=X[-1]+deltaX
        X=np.append(X,Xnew)
        #solve for Vz
        az=dVzdt(Ft[-1],Fd[-1],Fl[-1],M[-1],Vz[-1],Vx[-1],g[-1])
        az_array=np.append(az_array,az)
        Vznew=Vz[-1]+dt*az
        
        #Integrate using Trapezoidal rule 
        deltaZ=dt/2*(Vz[-1]+Vznew)
        Znew=Z[-1]+deltaZ
        Z=np.append(Z,Znew)
        print(Znew)
        #update velocities
        Vx=np.append(Vx,Vxnew)
        Vz=np.append(Vz,Vznew)
        
        
        #######################################################
        #   Update Parameters                                 #
        #######################################################
        V=np.append(V,np.sqrt(Vxnew*Vxnew+Vznew*Vznew))
        Machnew=get_Mach(V[-1],a) 
        Mach=np.append(Mach,Machnew)
        gnew=get_g(mu,Req,R,Z[-1],i0,J2)
        g=np.append(g,gnew)
        p=np.append(p,get_p(Z[-1]))
        T=np.append(T,get_T(Z[-1]))
        rho=np.append(rho,get_rho(p[-1],T[-1],Rgas))
        Clnew,Cdnew=H_aerodynamics_coefficients(Mach[-1],0)
        Cl=np.append(Cl,Clnew)
        Cd=np.append(Cd,Cdnew)
        Fdnew=get_Fd(Cdnew,rho[-1],S,V[-1])
        Fd=np.append(Fd,Fdnew)
        Flnew=get_Fl(Clnew,rho[-1],S,V[-1])
        Fl=np.append(Fl,Flnew)
        M=np.append(M, M[-1] - mdot[-1]*dt)
        
        if t_tot<=tb:
            #TWratio=np.append(TWratio,get_TWprofile(mode,t_tot,g[-1],1.5,4,tb,i0,Z[-1]))
            mdotnew=1/ceff*(M[-1]*g[-1]*TWratio[i]-Ae*(pe-p[-1]))
            Ftnew=get_Ft(mdot[-1],ceff,Ae,pe,p[-1])
            
        else:
            #TWratio=np.append(TWratio,0)
            mdotnew=0
            Ftnew=0
            
        Ft=np.append(Ft,Ftnew)
        mdot=np.append(mdot,mdotnew)


t_array=np.linspace(0,t_tot,len(Z)) 
a_array=np.sqrt(ax_array*ax_array+az_array*az_array) 
gamma=np.arctan(Vz/Vx)*180/np.pi                  #in degrees

print("~~~ Attained percentage of phasing velocity: ", (V[-1]+Vxfree)/V_phasing*100," % ~~~") #accounts for Rotatio of Mars
print("~~~ Attained final flight path angle: ",gamma[-1]," deg ~~~")


Mprop=np.sum(mdot*dt)
print("~~~ Propellant mass needed: ",Mprop," kg ~~~")
ascent_DeltaV=ceff*np.log(M[0]/(M[0]-Mprop))
print("~~~ DeltaV needed: ", ascent_DeltaV," m/s ~~~")

q=0.5*rho*V*V
#=====================================================================================================================================================================================================================
# Verification
#=====================================================================================================================================================================================================================
dV_verification=ROM_deltaV(g[0],R,h_phasing,h0,omega_mars,i0,V_phasing,tb)
#print(dV_verification)


#=====================================================================================================================================================================================================================
# Plotting
#=====================================================================================================================================================================================================================
if plotting:
    
    #time vs Z
    plt.plot(t_array,Z/10**3,color="navy")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Altitude")
    plt.xlabel("Time [s]")
    plt.ylabel("Altitude [km]")
    plt.show()
    
    #X vs Z
    plt.figure()
    plt.plot(X/10**3,Z/10**3,color="red")    
    plt.grid(color="gainsboro")
    plt.title("Lateral distance vs Altitude")
    plt.xlabel("Lateral distance [km]")
    plt.ylabel("Altitude [km]")
    plt.show()
    
    #velocity vs Z
    plt.figure()
    plt.plot(V/10**3,Z/10**3,color="orange")    
    plt.grid(color="gainsboro")
    plt.title("Velocity vs Altitude")
    plt.xlabel("Velocity [km/s]")
    plt.ylabel("Altitude [km]")
    plt.show()

    #Z vs gamma
    plt.figure()
    plt.plot(Z/10**3,gamma,color="hotpink")    
    plt.grid(color="gainsboro")
    plt.title("Altitude vs Flight Path angle")
    plt.xlabel("Altitude [km]")
    plt.ylabel("Flight Path angle [deg]")
    plt.show()

    #t vs Ft
    plt.figure()
    plt.plot(t_array,Ft/10**3,color="lime")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Thrust")
    plt.xlabel("Time [s]")
    plt.ylabel("Thrust [kN]")
    plt.show()

    #t vs mdot
    plt.figure()
    plt.plot(t_array,mdot,color="gold")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Mass Flow")
    plt.xlabel("Time [s]")
    plt.ylabel("Mass Flow [kg/s]")
    plt.show()

    #z vs q
    plt.figure()
    plt.plot(Z/10**3,q,color="cyan")    
    plt.grid(color="gainsboro")
    plt.title("Altitude vs Dynamic Pressure")
    plt.xlabel("Altitude [km]")
    plt.ylabel("Dynamic Pressure [Pa]")
    plt.show()









print("[---", (time.time() - start_time) ,"seconds ---]" )