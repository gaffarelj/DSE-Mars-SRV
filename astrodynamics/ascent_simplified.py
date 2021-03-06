# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:51:19 2020

@author: lucat
"""

import numpy as np
import matplotlib.pyplot as plt
from BAT import get_a, get_Mach, get_g, get_p, get_T, get_rho, get_ceff, Lambda, Mprop, get_Fl, get_Fd, get_Ft, get_TWratio, massfractions,ROM_deltaV,Mpvertical
#from Aero.aero_calcs import aerodynamics_coefficients
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
from aero_calcs import SS_aerodynamics_coefficients
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
gamma_gas=1.37                  #[-] heat capacity of Martian air
#=====================================================================================================================================================================================================================
# Node's properties & phasing orbit & Mars base
#=====================================================================================================================================================================================================================
i_node=42.48*np.pi/180             #[rad]
h_node=500*10**3                #[m]
i_phasing=i_node                #[rad]
V_phasing=3272.466              #[m/s]
h_phasing=609.74*10**3          #[m]
h0=-3*10**3                     #[m] altitude of the launch site wrt volumetric mean altitude
i0=i_node                       #[rad]
#=====================================================================================================================================================================================================================
# Vehicle's properties
#=====================================================================================================================================================================================================================
Isp=383.250565907662                       #[s] LOX-LCH4 383.250565907662
ceff=get_ceff(Isp)             #[m/s]

d=6                            #[m] diameter of the vehicle assuming cylinder
n=9                            #[-] number of engines
S=np.pi/4*d**2                 #[m^2] reference surface area used here is the area of the circle of the cross-section
De=1.35049466031671                          #[m] Diameter of an engine
Ae=np.pi/4*De*De               #[m] Exit area of engine
pe=6077.910186177842                     #[Pa] exit pressure of exhaust gasses
                      


#=====================================================================================================================================================================================================================
# Switches
#=====================================================================================================================================================================================================================

plotting=True      #Do you wanna plot? no=False
updateMOI=True
#=====================================================================================================================================================================================================================
# Simplified 2D point-mass model: gravity turn
#=====================================================================================================================================================================================================================
def TWparab(t,t0,tb,TW0,TWe):
	
	t_list=[t0,(tb+t0)/3,tb]
	TW_list=[TW0,2.85,TWe]
	
	def parabola(x, a, b, c):
		return a*x**2 + b*x + c
	
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_list, TW_list)
	TWparab=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	return TWparab

def TWlinear(t,t0,tb,TW0,TWe):

    TWlinear=(TWe-TW0)/tb*t+TW0

    return TWlinear


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

def get_Mp(DeltaV,ceff,Mwet):
	return (Mwet*(np.exp(DeltaV/ceff)-1))/(np.exp(DeltaV/ceff))

#initial conditions
#V0=8.33
V0=110/3.6                                               #8.33[m/s] circa 30 km/h
gamma0=(90-0.51)*np.pi/180                               #degs off vertical axis 3.2


Vx0=V0*np.cos(gamma0)
Vz0=V0*np.sin(gamma0)
#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
Mp, Mwet=massfractions(ceff,202413.4011,4495.16936199617+91.66666666666633,45.842,289.485,88,30,391.492,373.65,487.238)  #old Mwet:198948 kg, old DV_ascent: 4188.466 m/s
M=np.array([Mwet[0]])
print(M)

    
Vx=np.array([Vx0])
Vz=np.array([Vz0])
V=np.array([V0])
    
#earth reference frame is on the surface of Mars (where the base is, i.e. -3km)
#Z=np.array([10])    #Z0=10m
#X=np.array([0.524]) #X0=0.524m
Z=np.array([100])
X=np.array([0])

TW0=1.5
TWe=4		 #4

tb=489.5180709799191  #489.5180709799191 ceff/(TW0*9.80665)*np.log(M[0]/(M[0]-Mp[0]))+185

Mpver=Mpvertical(ceff,Mwet[0],TW0,V0)
DeltaVver=ceff*np.log(Mwet[0]/(Mwet[0]-Mpver))
tbver=ceff/(3.71*TW0)*np.log(Mwet[0]/(Mwet[0]-Mpver))
M=np.array([Mwet[0]-Mpver])
    
g=np.array([get_g(mu,Req,R,Z[-1],i0,J2)])
dt=0.01
t_tot=0
    
#*9.80665/g0_mars 

p=np.array([get_p(Z[-1])])
T=np.array([get_T(Z[-1])])
rho=np.array([get_rho(p[-1],T[-1],Rgas)])
Mach=np.array([get_Mach(V[0],get_a(gamma_gas,Rgas, T[0]))])
#TWratio=np.array([get_TWprofile(mode,0,g[0],1.5,4,tb,i0,Z[0])])
#TWratio=[TWparab(0,0,tb,TW0,TWe)]
TWratio=[TWlinear(0,0,tb,TW0,TWe)]
Cl0,Cd0=SS_aerodynamics_coefficients(Mach,0)
Cl=np.array([Cl0])
Cd=np.array([Cd0])
Fd=np.array([get_Fd(Cd[0],rho[0],S,V[0])])
Fl=np.array([get_Fl(Cl[0],rho[0],S,V[0])])
Ft=np.array([TWratio[0]*M[0]*g[0]])
mdot=np.array([(Ft[0]-Ae*(pe-p[0]))/(ceff)])
ax_array=np.array([dVxdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1])])
az_array=np.array([dVzdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1],g[-1])])
#az_array=np.array([0])
#ax_array=np.array([0])

#TWratio=np.array([Ft[0]/(M[0]*g[0])])
    
gamma=[gamma0]
gamma_dot=[0]    
while Z[-1]<h_phasing+abs(h0) and Z[-1]>=0:
#for i in range(2):
        
        #######################################################
        #   Compute parameters needed                         #
        #######################################################
        t_tot+=dt 
        
    
        
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
        gammanew=np.arctan(Vznew/Vxnew)
        gamma.append(gammanew)
        
        
        #######################################################
        #   Update Parameters                                 #
        #######################################################
        V=np.append(V,np.sqrt(Vxnew*Vxnew+Vznew*Vznew))
        gamma_dot.append(-g[-1]/V[-1]*np.cos(float(gammanew))+float(Fl[-1])/M[-1]*np.cos(float(gammanew)))

        p=np.append(p,get_p(Z[-1]))
        T=np.append(T,get_T(Z[-1]))
        a=get_a(gamma_gas,Rgas, T[-1])
        Machnew=get_Mach(V[-1],a) 
        Mach=np.append(Mach,Machnew)
        gnew=get_g(mu,Req,R,Z[-1],i0,J2)
        g=np.append(g,gnew)        
        rho=np.append(rho,get_rho(p[-1],T[-1],Rgas))
        Clnew,Cdnew=SS_aerodynamics_coefficients(Mach[-1],0)
        Cl=np.append(Cl,Clnew)
        Cd=np.append(Cd,Cdnew)
        Fdnew=get_Fd(Cdnew,rho[-1],S,V[-1])
        Fd=np.append(Fd,Fdnew)
        Flnew=get_Fl(Clnew,rho[-1],S,V[-1])
        Fl=np.append(Fl,Flnew)
        M=np.append(M, M[-1] - mdot[-1]*dt)
        
        if t_tot<=tb:
            #TWratio=np.append(TWratio,get_TWprofile(mode,t_tot,g[-1],1.5,4,tb,i0,Z[-1]))            
            """acc=np.sqrt(ax_array[-1]*ax_array[-1]+az_array[-1]*az_array[-1])
            if acc>=4*9.81/3.71:
               Ftnew=(0.002695*np.sqrt(-137641*((g[-1]**2*Vx[-1]**2-111.869*(V[-1]**2))*M[-1]**2*V[-1]**2-2*g[-1]*Fl[-1]*M[-1]*Vx[-1]*(V[-1]**2)*V[-1]+Fl[-1]**2*(V[-1]**2)**2))+371*(g[-1]*M[-1]*Vz[-1]*V[-1]+Fd[-1]*V[-1]**2))/V[-1]**2
               mdotnew=(Ftnew-Ae*(pe-p[-1]))/(ceff)
            """
            TWratio.append(TWlinear(t_tot,0,tb,TW0,TWe))
            Ftnew=TWratio[-1]*M[-1]*g[-1]
            mdotnew=(Ftnew-Ae*(pe-p[-1]))/(ceff)
        else:
            #TWratio=np.append(TWratio,0)
            TWratio.append(0)
            Ftnew=0
            mdotnew=0
            
            
        Ft=np.append(Ft,Ftnew)
        mdot=np.append(mdot,mdotnew)



t_array=np.linspace(0,t_tot,len(Z)) 
a_array=np.sqrt(ax_array*ax_array+az_array*az_array)    
gamma=np.array(gamma)*180/np.pi

gamma_dot=np.array(gamma_dot)
print("~~~ Attained percentage of phasing velocity: ", (V[-1]+Vxfree)/V_phasing*100," % ~~~") #accounts for Rotatio of Mars
print("~~~ Attained final flight path angle: ",gamma[-1]," deg ~~~")


Mprop=np.sum(mdot*dt)
print("~~~ Propellant mass needed: ",Mprop," kg ~~~")
ascent_DeltaV=ceff*np.log(M[0]/(M[0]-Mprop))
print("~~~ DeltaV needed: ", ascent_DeltaV," m/s ~~~")

q=0.5*rho*V*V

#Save the pitch rate gamma_dot from the gravity turn as a txt file
np.savetxt("gravity_pitchrate.txt", list(gamma_dot), delimiter=",")


ij=list(a_array).index(np.max(a_array))
a_crew=list(a_array[:ij]-1)+list(np.zeros((len(a_array)-ij,1)))
a_crew=np.array(a_crew)
#=====================================================================================================================================================================================================================
# Verification
#=====================================================================================================================================================================================================================
dV_verification=ROM_deltaV(g[0],R,h_phasing,h0,omega_mars,i0,V_phasing,tb)
#print(dV_verification)


#=====================================================================================================================================================================================================================
# Plotting
#=====================================================================================================================================================================================================================
if plotting:
	
	#time vs a_crew
    plt.plot(t_array,(a_array+3.71)/9.81,color="navy")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Acceleration experienced by the Crew")
    plt.xlabel("Time [s]")
    plt.ylabel("Acceleration [g_earth's]")
    plt.show()
	
	#a vs time
    plt.figure()
    plt.plot(t_array,a_array/9.81,color="cyan")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Acceleration")
    plt.xlabel("Time [s]")
    plt.ylabel("Acceleration [g_earth's]")
    plt.show()
"""
    #time vs Z
    plt.figure()
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
    plt.plot(t_array,Ft*1/10**3,color="lime")    
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


    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('Altitude [km]', color="navy")
    ax1.plot(t_array,Z/10**3,color="navy", label="Altitude")
    ax1.tick_params(axis='y', labelcolor="navy")
    plt.legend(loc="upper left")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Velocity [km/s]', color="darkorange")
    ax2.plot(t_array, V/10**3, color="darkorange", label="Velocity")
    ax2.tick_params(axis='y', labelcolor="darkorange")
    plt.grid(color="gainsboro")
    plt.title("Altitude and Velocity vs time")
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend(loc="upper right")
    plt.show()


    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('Acceleration [g_earth]', color="red")
    ax1.plot(t_array,a_crew*3.71/9.81,color="red", label="Acceleration")
    ax1.tick_params(axis='y', labelcolor="red")
    plt.grid(axis='y',color="gainsboro")
    plt.legend(loc="upper left")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Thrust [kN]', color="green")
    ax2.plot(t_array,Ft*1/10**3 , color="green", label="Thrust")
    ax2.tick_params(axis='y', labelcolor="green")
    
    plt.title("Acceleration and Thrust vs time")
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend(loc="upper right")
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('Flight path angle [deg]', color="hotpink")
    ax1.plot(t_array,gamma,color="hotpink", label="Flight path angle")
    ax1.tick_params(axis='y', labelcolor="hotpink")
    plt.grid(axis='y',color="gainsboro")
    plt.legend(loc="upper left")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Mass flow [kg/s]', color="dodgerblue")
    ax2.plot(t_array,mdot , color="dodgerblue", label="Mass flow")
    ax2.tick_params(axis='y', labelcolor="dodgerblue")
    
    plt.title("Flight path angle and Mass flow vs time")
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend(loc="upper right")
    plt.show()

print("[---", (time.time() - start_time) ,"seconds ---]" )




	
	
	
	
	
	
	
Mdry=33293.30437	
if updateMOI:
	Mpremain0=M[0]-Mdry
	i=list(mdot).index(np.max(mdot))
	Mp1=np.sum(mdot[:i]*dt)
	Mpremain1=M[0]-Mp1-Mdry
	#final mass flow point
	j=list(mdot).index(mdot[mdot!=0][-1])
	Mp2=np.sum(mdot[:j]*dt)
	Mpremain2=M[0]-Mp2-Mdry
	#point between max mass flow and final mass flow point
	lm=round((i+j)/2)
	Mp3=np.sum(mdot[:lm]*dt)
	MpremainP1=M[0]-Mp3-Mdry	
	#Point between start of gravity turn and mdot max. point:
	km=round((i)/2)
	Mp4=np.sum(mdot[:km]*dt)
	MpremainP2=M[0]-Mp4-Mdry
	
	
	print()
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("Point of start of gravity turn: ",Mpremain0)
	print("Point between start of gravity turn and mdot max. point: ", MpremainP2)
	print("Point of mdot max.: ",Mpremain1)
	print("Point between mdot max. and mdot final: ",MpremainP1)
	print("Point of mdot final: ",Mpremain2)
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
"""


Mateusz=False
if Mateusz:
	k=list(a_array).index(np.max(a_array))
	ang1=np.arctan(az_array[k]/ax_array[k])
	a_long=np.max(a_array)*np.cos((gamma[k]*np.pi/180-ang1))*3.71/9.81
	a_lat=np.max(a_array)*np.sin((gamma[k]*np.pi/180-ang1))*3.71/9.81
	jj=list(Ft).index(np.max(Ft))
	ang2=np.arctan(az_array[jj]/ax_array[jj])
	a_long2=a_array[jj]*np.cos((gamma[jj]*np.pi/180-ang2))
	a_lat2=a_array[jj]*np.sin((gamma[jj]*np.pi/180-ang2))
	
	print()
	print("FOR MATEUSZ:")
	print("Longitudinal acceleration: ", a_long2, " in Mars g's")
	print("Lateral acceleration: ", a_lat2, " in Mars g's")	
	"""
	print()
	print("FOR MATEUSZ:")
	print("Longitudinal acceleration: ", a_long, " in earth g's")
	print("Lateral acceleration: ", a_lat, " in earth g's")
	print("Remaining overall mass at this point: ",M[k])
	print()
	"""
	ang1=np.arctan(Vz[k]/Vx[k])
	V_long=V[k]*np.cos((gamma[k]*np.pi/180-ang1))
	V_lat=V[k]*np.sin((gamma[k]*np.pi/180-ang1))	
	print("FOR DIMITRIS:")
	print("Longitudinal acceleration: ", V_long, " in earth g's")
	print("Lateral acceleration: ", V_lat, " in earth g's")
	print("Remaining overall mass at this point: ",M[k])