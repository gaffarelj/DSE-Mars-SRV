# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:06:20 2020

@author: lucat
"""

import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
from aero_calcs import H_aerodynamics_coefficients
#========================================================================================================================================================================================================================
# Functions
#========================================================================================================================================================================================================================
#Mach
def get_Mach(V,a):
	Mach = V / a
	return Mach
#speed of sound
def get_a(gamma,Rgas,T):
	a = np.sqrt(gamma * Rgas * T)
	return a
#mu=gravitational parameter Mars, Req=equatorial radius Mars, R=mean volumetric
#radius Mars, h=altitude, theta=latitude [rad]
#gravitational acceleration
def get_g(mu,Req,R,h,theta,J2):
	P2 = 3 / 2 * np.sin(theta) ** 2 - 1 / 2
	g = mu / (R + h) ** 2 * (1 - 3 * J2 * (Req / (R + h)) ** 2 * P2)
	return g

#From mars_standard_atmosphere code [based on Viking 1 measurements]:

#pressure
def get_p(h):

	if h < 39000:
		p = 610.5 * (228.5 / (228.5 - 1.8 * h / 1000)) ** (19.435 / -1.8)
	elif h < 48000:
		p = 11.6025 * np.exp(-19.435 * (h / 1000 - 39) / 158.3)
	elif h < 55000:
		p = 3.84305 * (158.3 / (158.3 - 2.35 * (h / 1000 - 48))) ** (19.435 / -2.35)
	elif h < 66000:
		p = 1.55091 * (141.85 / (141.85 + 0.65 * (h / 1000 - 55))) ** (19.435 / 0.65)
	elif h < 75000:
		p = 0.356464 * (149 / (149 - 2.5 * (h / 1000 - 66))) ** (19.435 / -2.5)
	elif h < 84000:
		p = 0.099843 * (126.5 / (126.5 + 2.5 * (h / 1000 - 75))) ** (19.435 / 2.5)
	elif h < 95000:
		p = 0.0279653 * np.exp(-19.435 * (h / 1000 - 84) / 149)
	elif h < 105000:
		p = 0.00666032 * (149 / (149 - 1.4 * (h / 1000 - 95))) ** (19.435 / -1.4)
	elif h < 115897:
		p = 0.00169282 * (135 / (135 - 0.65 * (h / 1000 - 105))) ** (19.435 / -0.65)
	else:
		p = 0
	return p 

#Temperature
def get_T(h):

	if h < 39000:
		T = 228.5 - 1.8 * h / 1000
	elif h < 48000:
		T = 158.3
	elif h < 55000:
		T = 271.1 - 2.35 * h / 1000
	elif h < 66000:
		T = 106.1 + 0.65 * h / 1000
	elif h < 75000:
		T = 314 - 2.5 * h / 1000
	elif h < 84000:
		T = -61 + 2.5 * h / 1000
	elif h < 95000:
		T = 149
	elif h < 105000:
		T = 282 - 1.4 * h / 1000
	elif h < 115897:
		T = 203.25 - 0.65 * h / 1000
	else:
		e = (h / 1000 - 120) * (3389.51 + 120) / (3389.51 + h / 1000)
		T = 200 - 72.225 * np.exp(-0.0195 * e)
	return T

#density
def get_rho(p, T, Rgas):
	return p / (Rgas * T)

#equivalent jet velocity
def get_ceff(Isp):
	ceff = 9.80665 * Isp
	return ceff

#Mass ratio.  Note M0=Mi+Mp
def Lambda(M0,Me):
	Lambda = M0 / Me
	return Lambda


#propellant mass: rewritten Tsiolkovksy in terms of wet mass
def Mprop(ceff,Mwet,deltaV):
	Mprop = (Mwet * (np.exp(deltaV / ceff) - 1)) / np.exp(deltaV / ceff)
	return Mprop


#thrust force, where: mdot=mass flow,ceff=equivalent jet velocity,Ae=nozzle
#exit area,pe=nozzle exit pressure,p=atmospheric pressure at given altitude
def get_Ft(mdot,ceff,Ae,pe,p):
	Ft = mdot * ceff + Ae * (pe - p)
	return Ft

#drag force
def get_Fd(Cd,rho,S,V):
	Fd = 0.5 * rho * V * V * S * Cd
	return Fd

#lift force
def get_Fl(Cl,rho,S,V):
	Fl = 0.5 * rho * V * V * S * Cl
	return Fl

#specific thrust
def psi_sp(Ft,M,g):
	psi_sp = Ft / (M * g)
	return psi_sp

#COmputes Thrust to weight ratio of the vehicle at given conditions
def get_TWratio(Ft,M,g):
	TWratio = Ft / (M * g)
	return TWratio

"""
def get_TWprofile(mode,t,gmars,TW0,TWe,tb,i,Z):
	J2earth=1082.63*10**(-6)
	Rearth=6371*10**3
	muearth=3.986004418*(10**14)
	TWearth=(TWe-TW0)/tb*t+TW0
	gearth=9.80665*Rearth*Rearth/(Rearth+Z)**2
	#gearth=muearth*(1/(Rearth+Z)**2-J2earth*Rearth*Rearth*3/(Rearth+Z)**4*(-0.5+3/2*np.sin(i)**2))
	print("this is earth g ", gearth)
	if mode=="linear":
		TWearth=(TWe-TW0)/tb*t+TW0
	TWmars=TWearth*gearth/gmars
	return TWmars
"""

#"phasing"=orbital velocity of phasing orbit + gravity loss
#"Hohmann"=hohmann transfer from phasing orbit to LMO
#"inclination"=inclination change of 5 deg at LMO
#"reentry"=burn to insert into orbit needed for reentry (apoasis 300km under R)
#"reserve"=10% reserves of tot deltav
def massfractions(ceff,M_wet,DV_ascent,DV_hohmann,DV_inclination,DV_rendezvous,DV_docking,DV_reentry,DV_landing,DV_reserve):
	#phase=["phasing","Hohmann","inclination","rendezvous","docking","reentry","landing","reserve"]
	#deltaV_c2=np.array([3272.466+916.29,
	#45.842,289.485,88,30,391.492,373.65,487.238])
	deltaV = np.array([DV_ascent,DV_hohmann,DV_inclination,DV_rendezvous,DV_docking,DV_reentry,DV_landing,DV_reserve])
	Mp = []
	Mwet = [M_wet]
	for i in range(5):
		
		Mp.append(Mprop(ceff,Mwet[i],deltaV[i]))
		Mwet.append(Mwet[i] - Mp[i])
		
	Mp = np.array(Mp)
	Mwet = np.array(Mwet)
	return Mp, Mwet


#Rough order of magnitude for ascent delta V
#h0=altitude of launch site, horbit=altitude of orbit to be attained,
#omega_mars=angular velocity of Mars, i0=inclination of the launch site
def ROM_deltaV(g0,R,h,h0,omega_mars,i0,Vorbit,tb):
	Vrot = omega_mars * (R + h0) * np.cos(i0)
	Vvert = np.sqrt(2 * g0 * h * (R) ** 2 / (R + h) ** 2)
	Vgloss = g0 * tb
	dV = np.sqrt((Vorbit - Vrot) ** 2 + Vvert ** 2) + Vgloss
	return dV

#acceleration in x
def dVxdt(T,D,L,M,Vz,Vx):
	dVxdt = (T - D) / M * Vx / (np.sqrt(Vx * Vx + Vz * Vz)) - L / M * Vz / (np.sqrt(Vx * Vx + Vz * Vz))
	return dVxdt

#acceleration in z
def dVzdt(T,D,L,M,Vz,Vx,g):
	dVzdt = (T - D) / M * Vz / (np.sqrt(Vx * Vx + Vz * Vz)) - g + L / M * Vx / (np.sqrt(Vx * Vx + Vz * Vz))
	return dVzdt

#point mass ascent simulation
def ascent_sim(tb,initial_tilt,i_base,h0,d,M_initial,Mp_class2,Isp,n,De,pe):
	"""
	Function that simulates the gravity turn ascent profile constrained by a linear variation of the T/W_ratio (T/W_0=1.5, T/W_final=4).
	It simulates it for a point mass in 2D accounting for aerodynamic forces. 
	Reference frame is on the surface of Mars, at the launch pad. Z is the altitude and X the lateral distance (points in the diraction of the orbit groundtrack).
	
	
	INPUTS:
		
		-tb: burn time [s]
	   
		-initial_tilt: the tilt in the initial flight path angle (at the launch pad) off the vertical axis [deg]. An initial_tilt=0
						would mean the rocket starts vertically.

		-i_base: latitude of the launch pad [deg]
		
		-h0: altitude of the launch pad wrt the volumetric mean altitude of Mars [m].
		
		-d: diameter of the vehicle [m]. Symmetric vehicle assumed
		
		-M_initial: mass of the vehicle sitting on the launch pad [kg]
		
		-Mp_class2: propellant Mass needed to attain phasing orbit, as estimated in the Class II deltaV for single stage vehicle [m/s].
		
		
		-Isp: specific impulse of the engines [s].
		
		-n: number of engines [-].
		
		-De: diameter of the exit plane of the nozzle [m].
		
		-pe: pressure of the exhaust gasses at the exti plane of the nozzle [Pa]
		
		-
	OUTPUTS:
		
		-V: array of the velocity throughout flight
		
		-Vxfree: float of the free velocity of the vehcile thanks to Mars rotation at that latitude
		
		-ascent_DeltaV: float of the predicted delta_V required to obtain the final conditions
	
		-q: array containing the dynamic pressure throughout flight
		
		-others: a dictionary containing other parameters that other people might need. 
		
		!NOTE: if there's another parameter that you need that was computed, please add it inside the "others" dictionary to return it.
		
	"""
	
	#=====================================================================================================================================================================================================================
	# Mars properties = Atmospheric composition: 95.32% CO2, 2.7% N2, 1.6% Ar,
	# 0.13% O2, 0.08% CO
	#=====================================================================================================================================================================================================================
	Req = 3396.2 * 10 ** 3                #[m] equatorial radius
	R = 3389.5 * 10 ** 3                  #[m] volumetric mean radius
	g0_mars = 3.71                    #[m/s^2] surface gravity
	mu = 0.042828 * 10 ** 6 * (10 ** 3) ** 3    #[m^3/s^2] gravitational parameter
	J2 = 1960.45 * 10 ** (-6)
	omega_mars = 0.7088218 * 10 ** (-4)   #[rad/s] Martian angular velocity
	M_atm = 43.34                     #[g/mol] Mean molecular weight
	Rgas = 8314.4621 / M_atm            #[J/(kg*K)] Martian air gas constant
	gamma_gas = 1.37                      #[-] heat capacity of Martian air
	#=====================================================================================================================================================================================================================
	# Node's properties & phasing orbit & Mars base
	#=====================================================================================================================================================================================================================
	i_node = i_base * np.pi / 180           #[rad]
	i_phasing = i_node                #[rad]
	V_phasing = 3272.466              #[m/s]
	h_phasing = 609.74 * 10 ** 3          #[m]
	#h0 [m] = altitude of the launch site wrt volumetric mean altitude.  -3km
	i0 = i_node                       #[rad]
	
	S = np.pi / 4 * d ** 2                 #[m^2] reference surface area used here is the area of the
										   #circle of the cross-section
	ceff = get_ceff(Isp)             #[m/s]
	Ae = np.pi / 4 * De * De               #[m] Exit area of engine
	
	#free Vx due to the rotation of Mars: Vx(t=0):
	Vxfree = omega_mars * (R + h0) * np.cos(i0)
	
	#initial conditions
	V0 = 0.00001                      #8.33[m/s] circa 30 km/h
	gamma0 = (90 - initial_tilt) * np.pi / 180      #degs off vertical axis 2.45
	
	Vx0 = V0 * np.cos(gamma0)
	Vz0 = V0 * np.sin(gamma0)
	#=====================================================================================================================================================================================================================
	# Simulation
	#=====================================================================================================================================================================================================================
	M = np.array([M_initial])
		
	Vx = np.array([Vx0])
	Vz = np.array([Vz0])
	V = np.array([V0])
		
	#earth reference frame is on the surface of Mars (where the base is, i.e.
	#-3km)
	Z = np.array([0.0])    #Z0=10m
	X = np.array([0.0]) #X0=0.524m
	
	TW0 = 1.5
	TWe = 4
		
	g = np.array([get_g(mu,Req,R,Z[-1],i0,J2)])
	dt = 0.01
	t_tot = 0
		
	 
	p = np.array([get_p(Z[-1])])
	T = np.array([get_T(Z[-1])])
	rho = np.array([get_rho(p[-1],T[-1],Rgas)])
	Mach = np.array([get_Mach(V[0],get_a(gamma_gas,Rgas, T[0]))])
	#TWratio=np.array([get_TWprofile(mode,0,g[0],1.5,4,tb,i0,Z[0])])
	TWratio = np.linspace(TW0 * 9.80665 / g0_mars,TWe * 9.80665 / g0_mars,round(tb / dt))
	Cl0,Cd0 = H_aerodynamics_coefficients(Mach,0)
	Cl = np.array([Cl0])
	Cd = np.array([Cd0])
	Fd = np.array([get_Fd(Cd[0],rho[0],S,V[0])])
	Fl = np.array([get_Fl(Cl[0],rho[0],S,V[0])])
	mdot = np.array([1 / ceff * (M[0] * g[0] * TWratio[0] - Ae * (pe - p[0]))])
	Ft = np.array([get_Ft(mdot[0],ceff,Ae,pe,p)])
	ax_array = np.array([dVxdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1])])
	az_array = np.array([dVzdt(Ft,get_Fd(Cd,rho,S,V[-1]),get_Fl(Cl,rho,S,V[-1]),M[-1],Vz[-1],Vx[-1],g[-1])])
	
	#TWratio=np.array([Ft[0]/(M[0]*g[0])])
		
	i = -1
		
	while Z[-1] < h_phasing + abs(h0) and Z[-1] >= 0:
	#for i in range(1):
			
			#######################################################
			#   Compute parameters needed #
			#######################################################
			t_tot+=dt 
			i+=1
			
			a = get_a(gamma_gas,Rgas, T[-1]) 
			#######################################################
			#   Solve for ax, az, Vx, Vz, X, Z #
			#######################################################
			#solve for Vx
			ax = dVxdt(Ft[-1],Fd[-1],Fl[-1],M[-1],Vz[-1],Vx[-1])
			ax_array = np.append(ax_array,ax)
			Vxnew = Vx[-1] + dt * ax
			
			#Integrate using Trapezoidal rule
			deltaX = dt / 2 * (Vx[-1] + Vxnew)
			Xnew = X[-1] + deltaX
			X = np.append(X,Xnew)
			#solve for Vz
			az = dVzdt(Ft[-1],Fd[-1],Fl[-1],M[-1],Vz[-1],Vx[-1],g[-1])
			az_array = np.append(az_array,az)
			Vznew = Vz[-1] + dt * az
			
			#Integrate using Trapezoidal rule
			deltaZ = dt / 2 * (Vz[-1] + Vznew)
			Znew = Z[-1] + deltaZ
			Z = np.append(Z,Znew)
			progress = round(Znew/h_phasing*100, 2)
			print(progress)
			#update velocities
			Vx = np.append(Vx,Vxnew)
			Vz = np.append(Vz,Vznew)
			
			
			#######################################################
			#   Update Parameters #
			#######################################################
			V = np.append(V,np.sqrt(Vxnew * Vxnew + Vznew * Vznew))
			Machnew = get_Mach(V[-1],a) 
			Mach = np.append(Mach,Machnew)
			gnew = get_g(mu,Req,R,Z[-1],i0,J2)
			g = np.append(g,gnew)
			p = np.append(p,get_p(Z[-1]))
			T = np.append(T,get_T(Z[-1]))
			rho = np.append(rho,get_rho(p[-1],T[-1],Rgas))
			Clnew,Cdnew = H_aerodynamics_coefficients(Mach[-1],0)
			Cl = np.append(Cl,Clnew)
			Cd = np.append(Cd,Cdnew)
			Fdnew = get_Fd(Cdnew,rho[-1],S,V[-1])
			Fd = np.append(Fd,Fdnew)
			Flnew = get_Fl(Clnew,rho[-1],S,V[-1])
			Fl = np.append(Fl,Flnew)
			M = np.append(M, M[-1] - mdot[-1] * dt)
			
			if t_tot <= tb:
				#TWratio=np.append(TWratio,get_TWprofile(mode,t_tot,g[-1],1.5,4,tb,i0,Z[-1]))
				mdotnew = 1 / ceff * (M[-1] * g[-1] * TWratio[i] - Ae * (pe - p[-1]))
				Ftnew = get_Ft(mdot[-1],ceff,Ae,pe,p[-1])
				
			else:
				#TWratio=np.append(TWratio,0)
				mdotnew = 0
				Ftnew = 0
				
			Ft = np.append(Ft,Ftnew)
			mdot = np.append(mdot,mdotnew)
	
	
	t_array = np.linspace(0,t_tot,len(Z)) 
	a_array = np.sqrt(ax_array * ax_array + az_array * az_array) 
	gamma = np.arctan(Vz / Vx) * 180 / np.pi                  #in degrees
	Mprop = np.sum(mdot * dt)
	ascent_DeltaV = ceff * np.log(M[0] / (M[0] - Mprop))
	q = 0.5 * rho * V * V
	
	others = {"Mach": Mach,
			"Lift": Fl,
			"Drag": Fd,
			"Thrust": Ft,
			"Mass": M,
			"V": V,
			"acceleration x": ax_array,
			"acceleration z": az_array,
			"altitude": Z,
			"density": rho}
	
	return V, Vxfree, ascent_DeltaV, q, others