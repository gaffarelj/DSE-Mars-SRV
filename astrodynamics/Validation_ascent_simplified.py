# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:41:31 2020

@author: lucat
"""

"""
Created on Fri May 15 20:51:19 2020

@author: lucat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from BAT import get_a, get_Mach, get_g, get_p, get_T, get_rho, get_ceff, Lambda, Mprop, get_Fl, get_Fd, get_Ft, get_TWratio, massfractions,ROM_deltaV,Mpvertical
#from Aero.aero_calcs import aerodynamics_coefficients
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
from aero_calcs import SS_aerodynamics_coefficients  
from ISA_DEF import ISA
import time
start_time = time.time()

#=====================================================================================================================================================================================================================
# Earth properties   
#=====================================================================================================================================================================================================================
Req=6378.137*10**3                #[m] equatorial radius
R=6371.000*10**3                  #[m] volumetric mean radius
g0_mars=9.80665                    #[m/s^2] surface gravity
mu=0.39860*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter  
J2=1082.63*10**(-6) 
i_mars=25.19*np.pi/180          #[rad] inlination of Mars' equator
omega_mars=7.2921159 *10**(-5)   #[rad/s] Martian angular velocity
M_atm=28.97                      #[g/mol] Mean molecular weight
Rgas=8314.4621/M_atm            #[J/(kg*K)] Martian air gas constant
gamma_gas=1.4                  #[-] heat capacity of Martian air
#=====================================================================================================================================================================================================================
# Node's properties & phasing orbit & Mars base
#=====================================================================================================================================================================================================================
i_node=51*np.pi/180             #[rad] 28.52*np.pi/180 
i_phasing=i_node                #[rad]
h_phasing=44726          #[m]
h0=11617                    #[m] altitude of the launch site wrt volumetric mean altitude
i0=i_node                       #[rad]
#=====================================================================================================================================================================================================================
# Vehicle's properties
#=====================================================================================================================================================================================================================



d=8.7                            #[m] diameter of the vehicle assuming cylinder   
S=np.pi/4*d**2                 #[m^2] reference surface area used here is the area of the circle of the cross-section
                  #[Pa] exit pressure of exhaust gasses
  
#SRB-> delivers 71.4% of sea level thrust                    
Isp_b=242.1	#was computed!
pe_b=0
Ae_b=0
tb_b=124
ceff_b=Isp_b*9.80665
#SSME
Isp_me=395.9	#was computed!
pe_me=27357.23772
Ae_me=77.45*0.053
ceff_me=Isp_me*9.80665
#=====================================================================================================================================================================================================================
# Switches
#=====================================================================================================================================================================================================================

plotting=False       #Do you wanna plot? no=False
updateMOI=False
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


def dVzdt(T,D,L,M,Vz,Vx,g):
    dVzdt=(T-D)/M*Vz/(np.sqrt(Vx*Vx+Vz*Vz))-g+L/M*Vx/(np.sqrt(Vx*Vx+Vz*Vz))
    return dVzdt



#initial conditions
#V0=8.33
V0=447.9341                                             #8.33[m/s] circa 30 km/h
gamma0=(90-33)*np.pi/180                               #degs off vertical axis 3.2


Vx0=V0*np.cos(gamma0)
Vz0=V0*np.sin(gamma0)
#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
#Mp, Mwet=massfractions(ceff,202413.4011,4495.16936199617+91.66666666666633,45.842,289.485,88,30,391.492,373.65,487.238)  #old Mwet:198948 kg, old DV_ascent: 4188.466 m/s
M=np.array([1376301])
print(M)

    
Vx=np.array([Vx0])
Vz=np.array([Vz0])
V=np.array([V0])
    
#earth reference frame is on the surface of Mars (where the base is, i.e. -3km)
#Z=np.array([10])    #Z0=10m
#X=np.array([0.524]) #X0=0.524m
Z=np.array([11617])
X=np.array([0])

TW0=2 #2 g's at this altitude
TWe=3

tb=tb_b-60
 
g=np.array([get_g(mu,Req,R,Z[-1],i0,J2)])
dt=0.001
t_tot=0
    
#*9.80665/g0_mars 
T0, p0 , rh0 = ISA(288.15,1.225,g[-1],Rgas,-0.0065,0,0.001,0.0028,0,-0.0028,-0.002,101325,Z[-1])
p=np.array([p0])
T=np.array([T0])
rho=np.array([rh0])
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
mdot_me=np.array([(0.286*Ft[0]-Ae_me*(pe_me-p[0]))/(ceff_me)])
mdot_b=np.array([((1-0.286)*Ft[0])/(ceff_b)])
ax_array=np.array([dVxdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1])])
az_array=np.array([dVzdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1],g[-1])])
#az_array=np.array([0])
#ax_array=np.array([0])

#TWratio=np.array([Ft[0]/(M[0]*g[0])])
    
j=[]
k=0    
while Z[-1]<h_phasing and Z[-1]>=0:
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
        
        
        #######################################################
        #   Update Parameters                                 #
        #######################################################
        V=np.append(V,np.sqrt(Vxnew*Vxnew+Vznew*Vznew))
        Tnew, pnew , rhnew = ISA(288.15,1.225,g[-1],Rgas,-0.0065,0,0.001,0.0028,0,-0.0028,-0.002,101325,Z[-1])
        p=np.append(p,pnew)
        T=np.append(T,Tnew)
        a=get_a(gamma_gas,Rgas, T[-1])
        Machnew=get_Mach(V[-1],a) 
        Mach=np.append(Mach,Machnew)
        gnew=get_g(mu,Req,R,Z[-1],i0,J2)
        g=np.append(g,gnew)        
        rho=np.append(rho,rhnew)
        Clnew,Cdnew=SS_aerodynamics_coefficients(Mach[-1],0)
        Cl=np.append(Cl,Clnew)
        Cd=np.append(Cd,Cdnew)
        Fdnew=get_Fd(Cdnew,rho[-1],S,V[-1])
        Fd=np.append(Fd,Fdnew)
        Flnew=get_Fl(Clnew,rho[-1],S,V[-1])
        Fl=np.append(Fl,Flnew)
        M=np.append(M, M[-1] - mdot_me[-1]*dt-mdot_b[-1]*dt)

        if t_tot<=tb:
            #TWratio=np.append(TWratio,get_TWprofile(mode,t_tot,g[-1],1.5,4,tb,i0,Z[-1]))            
            """acc=np.sqrt(ax_array[-1]*ax_array[-1]+az_array[-1]*az_array[-1])
            if acc>=4*9.81/3.71:
               Ftnew=(0.002695*np.sqrt(-137641*((g[-1]**2*Vx[-1]**2-111.869*(V[-1]**2))*M[-1]**2*V[-1]**2-2*g[-1]*Fl[-1]*M[-1]*Vx[-1]*(V[-1]**2)*V[-1]+Fl[-1]**2*(V[-1]**2)**2))+371*(g[-1]*M[-1]*Vz[-1]*V[-1]+Fd[-1]*V[-1]**2))/V[-1]**2
               mdotnew=(Ftnew-Ae*(pe-p[-1]))/(ceff)
            """
            TWratio.append(TWlinear(t_tot,0,tb,TW0,TWe))
            Ftnew=TWratio[-1]*M[-1]*g[-1]
            mdotnew_b=(0.714*Ftnew)/(ceff_b)
            mdotnew_me=(0.286*Ftnew-Ae_me*(pe_me-p[-1]))/(ceff_me)
        else:
            #TWratio=np.append(TWratio,0)
            TWratio.append(0)
            Ftnew=0
            mdotnew_b=0
            mdotnew_me=0
        k+=1    
        if t_tot>10 and t_tot<10+2*dt:
            j.append(["10",k])
        elif t_tot>20 and t_tot<20+2*dt:
            j.append(["20",k])
        elif t_tot>30 and t_tot<30+2*dt:
            j.append(["30",k])			
        elif t_tot>40 and t_tot<40+2*dt:
            j.append(["40",k])		
        elif t_tot>50 and t_tot<50+2*dt:
            j.append(["50",k])	
        elif t_tot>60 and t_tot<60+2*dt:
            j.append(["60",k])			
			
        Ft=np.append(Ft,Ftnew)
        mdot_me=np.append(mdot_me,mdotnew_me)
        mdot_b=np.append(mdot_b,mdotnew_b)



t_array=np.linspace(0,t_tot,len(Z)) 
a_array=np.sqrt(ax_array*ax_array+az_array*az_array) 
gamma=np.arctan(Vz/Vx)*180/np.pi                  #in degrees




q=0.5*rho*V*V



#=====================================================================================================================================================================================================================
# Validation: using STS-121 Ascent Data
#=====================================================================================================================================================================================================================
#start at altitude after 10km and until before 50 km
t_STS=np.array([0,10,20,30,40,50,60])
M_STS=np.array([1376301,1277921,1177704,1075683,991872,913254,880377])
h_STS=np.array([11617,15380,19872,25608,31412,38309,44726])

def hSTSfit(t,t_STS,h_STS):	
	def parabola(x, a, b, c):
		   return a*x**2 + b*x + c
	   
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_STS, h_STS)
	h_STS_fit=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	perr = np.sqrt(np.diag(pcov))
	h_STS_fit_down=parabola(t,fit_params[0]-perr[0],fit_params[1]-perr[1],fit_params[2]-perr[2])
	h_STS_fit_up=parabola(t,fit_params[0]+perr[0],fit_params[1]+perr[1],fit_params[2]+perr[2])
	return h_STS_fit, h_STS_fit_down, h_STS_fit_up

def MSTSfit(t,t_STS,M_STS):	
	def parabola(x, a, b, c):
		   return a*x**2 + b*x + c
	   
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_STS, M_STS)
	M_STS_fit=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	perr = np.sqrt(np.diag(pcov))
	M_STS_fit_down=parabola(t,fit_params[0]-perr[0],fit_params[1]-perr[1],fit_params[2]-perr[2])
	M_STS_fit_up=parabola(t,fit_params[0]+perr[0],fit_params[1]+perr[1],fit_params[2]+perr[2])
	return M_STS_fit, M_STS_fit_down, M_STS_fit_up


def MAPE(At,Ft):
	""" Computes the Mean Absolute Percentage Error (MAPE), for a list of estimated values (Ft), wrt actual values (At)
	"""
	#actual value
	#number of points
	n=1/len(At)
	#estimated value
	diff=[]
	for i in range(len(At)):
		diff.append(np.abs(At[i]-Ft[i])/At[i])
		
	MAPE=100*n*np.sum(diff)
	return MAPE

Ft_Z=[Z[0],Z[10001],Z[20000],Z[30000],Z[40001],Z[50001],Z[60001]]
Ft_M=[M[0],M[10001],M[20000],M[30000],M[40001],M[50001],M[60001]]

Z_MAPE=MAPE(h_STS,Ft_Z)
M_MAPE=MAPE(M_STS,Ft_M)

h_STS_fit, h_STS_fit_down, h_STS_fit_up=hSTSfit(t_STS,t_STS,h_STS)
M_STS_fit, M_STS_fit_down, M_STS_fit_up=MSTSfit(t_STS,t_STS,M_STS)

#plt.scatter(t_STS,h_STS,color="red",marker="+")
plt.plot(t_STS,h_STS_fit,color="red",label="STS-121 parabolic fit")
plt.scatter(t_STS,h_STS,color="red",marker="+")
plt.fill_between(t_STS, h_STS_fit_down,h_STS_fit_up, color="gray", alpha=0.2)

#time vs Z

plt.plot(t_array,Z,color="navy",label="simulation")
plt.grid(color="gainsboro")
plt.title("Time vs Altitude: MAPE=" + str(round(Z_MAPE,2)) + "%")
plt.xlabel("Time [s]")
plt.ylabel("Altitude [m]")
plt.legend()
plt.show()

#plt.scatter(t_STS,h_STS,color="red",marker="+")
plt.figure()
plt.plot(t_STS,M_STS_fit,color="forestgreen",label="STS-121 parabolic fit")
plt.scatter(t_STS,M_STS,color="forestgreen",marker="+")
plt.fill_between(t_STS, M_STS_fit_down,M_STS_fit_up, color="gray", alpha=0.2)

#time vs Z

plt.plot(t_array,M,color="gold",label="simulation")
plt.grid(color="gainsboro")
plt.title("Time vs Mass: MAPE="+ str(round(M_MAPE,2)) + "%")
plt.xlabel("Time [s]")
plt.ylabel("Mass [kg]")
plt.legend()
plt.show()
        

	

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

    #a vs time
    plt.figure()
    plt.plot(t_array,a_array*3.71/9.81,color="cyan")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Acceleration")
    plt.xlabel("Time [s]")
    plt.ylabel("Acceleration [g_earth's]")
    plt.show()







print("[---", (time.time() - start_time) ,"seconds ---]" )


