# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 16:27:51 2020

@author: lucat
"""

import numpy as np
import scipy.optimize
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from BAT import get_a, get_Mach, get_g, get_p, get_T, get_rho, get_ceff, Lambda, Mprop, get_Fl, get_Fd, get_Ft, get_TWratio, massfractions,ROM_deltaV, ascent_sim
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
from aero_calcs import SS_aerodynamics_coefficients
import time
start_time = time.time()
#========================================================================================================================================================================================================================================================
#   Switches
#========================================================================================================================================================================================================================================================
plotting=True             #no plots? -> plotting=False
updateMOI=False           #need updated MOI and x_c.g. location?->updateMOI=True & then give to Dimitris

#========================================================================================================================================================================================================================================================
#   Functions 
#========================================================================================================================================================================================================================================================
def inputDimitris(Mprop,Mptot,mdot):
	"""Function that returns the two points (total propellant mass) for which I need the MOIs and c.g. xc.g. location
	   These points, in addition to the wet mass one, are: at max. mass flow point and end of burn time point"
	"""	
	#max.mass flow point
	
	i=list(mdot).index(np.max(mdot))
	Mp1=Mprop[i]
	Mpremain1=Mptot-(Mprop[0]-Mprop[i])
	#final mass flow point
	j=list(mdot).index(mdot[mdot!=0][-1])
	Mp2=Mprop[j]
	Mpremain2=Mptot-(Mprop[0]-Mprop[j])
	#random point before mdot_max
	k=round((i+0)/2)
	MpremainP1=Mptot-(Mprop[0]-Mprop[k])
	#random point after mdot_max
	l=round((i+j)/2)
	MpremainP2=Mptot-(Mprop[0]-Mprop[l])
	
	return Mpremain1, Mpremain2, MpremainP1, MpremainP2

def xcg(t,xcg0,t_xcg0,xcgm,t_xcgm,xcge,t_xcge,tb):
	"""computes the c.g. shift for a given time. In the reference frame of "vehicle layout"
		Fitted by parabola
	"""
	t_list=[t_xcg0,t_xcgm,t_xcge]
	xcg_list=[xcg0,xcgm,xcge]
	def parabola(x, a, b, c):
		return a*x**2 + b*x + c
	def ln(x,a,b):
		return a+b*np.log(x)
	def ex(x,a,b,c):
		return a*np.exp(b*x)+c 
	# calculate polynomial, third element is the order of the polynomial
	#z = np.polyfit(t_list, xcg_list, 2)
	#f = np.poly1d(z)
	#xcg_current=f(t)
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_list, xcg_list)
	xcg_current=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	
	return xcg_current
	
def delta_cg(t,xcg0,t_xcg0,xcgm,t_xcgm,xcge,t_xcge,tb):
	""" Computes the cg location through-out flight. 
	    Parabolic variation
	"""
	t_list=[t_xcg0,t_xcgm,t_xcge]
	xcg_list=[xcg0,xcgm,xcge]
	
	def parabola(x, a, b, c):
		return a*x**2 + b*x + c
	def ln(x,a,b):
		return a+b*np.log(x)
	def ex(x,a,b,c):
		return a*np.exp(b*x)+c 
	
	# calculate polynomial, third element is the order of the polynomial
	#z = np.polyfit(t_list, xcg_list, 2)
	#f = np.poly1d(z)
	#xcg_current=f(t)
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_list, xcg_list)
	xcg_current=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	delta_cg=xcg0-xcg_current
	#delta_cg=xcg0-xcg_current
	return delta_cg

def get_xcp(lcap,lcyl,ltot,xcg0,Acap,Acyl,delta_cg):
	"""Computes x-location of the center of pressure using the Projected Area Method
		
	"""
	xcp0=(Acap*1/4*lcap+Acyl*(lcap+0.5*lcyl))/(Acap+Acyl)
	xcp=ltot-xcg0-xcp0+delta_cg
	return xcp


def get_Ix(t,Ix0,t_xcg0,Ixm,t_xcgm,Ixe,t_xcge,tb):
	""" Computes mass moment of inertia in x direction.parabolic fit & derivative"""
	#Ix=(Ixe-Ix0)/tb*t+Ix0
	t_list=[t_xcg0,t_xcgm,t_xcge]
	Ix_list=[Ix0,Ixm,Ixe]
	
	def parabola(x, a, b, c):
		return a*x**2 + b*x + c
	def ln(x,a,b):
		return a+b*np.log(x)
	def ex(x,a,b,c):
		return a*np.exp(b*x)+c 
	
	# calculate polynomial, third element is the order of the polynomial
	#z = np.polyfit(t_list, Ix_list, 2)
	#f = np.poly1d(z)
	#Ix=f(t)
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_list, Ix_list)
	Ix=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	dIxdt=2*fit_params[0]*t+fit_params[1]
	#dIxdt=(f(t+0.1)-Ix)/0.1
	return [Ix,dIxdt]


def get_Iy(t,Iy0,t_xcg0,Iym,t_xcgm,Iye,t_xcge,tb):
	""" Computes mass moment of inertia in y direction & derivative"""
	#Ix=(Iye-Iy0)/tb*t+Iy0
	t_list=[t_xcg0,t_xcgm,t_xcge]
	Iy_list=[Iy0,Iym,Iye]
	
	def parabola(x, a, b, c):
		return a*x**2 + b*x + c
	def ln(x,a,b):
		return a+b*np.log(x)
	def ex(x,a,b,c):
		return a*np.exp(b*x)+c 
	
	# calculate polynomial, third element is the order of the polynomial
	#z = np.polyfit(t_list, Iy_list, 2)
	#f = np.poly1d(z)
	#Iy=f(t)
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_list, Iy_list)
	Iy=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	#Iy=f(t)
	dIydt=2*fit_params[0]*t+fit_params[1]
	#dIydt=(f(t+0.1)-Iy)/0.1
	return [Iy,dIydt]

def get_Iz(t,Iz0,t_xcg0,Izm,t_xcgm,Ize,t_xcge,tb):
	""" Computes mass moment of inertia in z direction & derivative"""
	#Iz=(Ize-Iz0)/tb*t+Iz0
	t_list=[t_xcg0,t_xcgm,t_xcge]
	Iz_list=[Iz0,Izm,Ize]
	
	def parabola(x, a, b, c):
		return a*x**2 + b*x + c
	def ln(x,a,b):
		return a+b*np.log(x)
	def ex(x,a,b,c):
		return a*np.exp(b*x)+c 
	
	# calculate polynomial, third element is the order of the polynomial
	#z = np.polyfit(t_list, Iz_list,2)
	#f = np.poly1d(z)
	#Iz=f(t)
	fit_params, pcov = scipy.optimize.curve_fit(parabola, t_list, Iz_list)
	Iz=parabola(t,fit_params[0],fit_params[1],fit_params[2])
	
	dIzdt=2*fit_params[0]*t+fit_params[1]
	#dIzdt=(f(t+0.1)-Iz)/0.1
	return [Iz,dIzdt]


def TWlinear(t,t0,tb,TW0,TWe):

    TWlinear=(TWe-TW0)/tb*t+TW0

    return TWlinear

def T_Cb(delta,tau,theta,phi,psi):
	""" Transformation Matrix to go from Body Fixed (Fb) to inertial (FI) reference frame
	    delta=Latitude, tau=Longitude
		phi=Roll Angle, theta=Pitch Angle, psi=Yaw Angle
	"""	
	#Rotation needed to go from Inertial Frame to Mars Centered Frame
	#Omega_t=rotational speed of the planet about itself, t0=stellar time point (here used for the J2000 inertial frame)
#	T_CI=np.matrix([[ np.cos(Omega_t)*t0       ,     np.sin(Omega_t)*t0       ,     0            ],
#				    [-np.sin(Omega_t)*t0       ,     np.cos(Omega_t)*t0       ,     0            ],
#					[ 0                        ,     0                        ,     1            ]])
	
	#Rotation needed to go from Mars Centered to Vehicle Carried Normal Frame
	T_EC=np.matrix([[-np.sin(delta)*np.cos(tau),	-np.sin(delta)*np.sin(tau),		np.cos(delta)],
				    [-np.sin(tau)              ,	 np.cos(tau)              ,     0            ],
					[-np.cos(delta)*np.cos(tau),    -np.cos(delta)*np.sin(tau),    -np.sin(delta)]])
	
	#Rotations to go from Vehicle Carried Normal frame to Body frame
	T_X=np.matrix([[  1 					   ,     0                        ,     0            ],
				   [  0                        ,     np.cos(phi)	          ,     np.sin(phi)	 ],
				   [  0                        ,    -np.sin(phi)              ,     np.cos(phi)  ]])
	
	T_Y=np.matrix([[  np.cos(theta)            ,     0                        ,    -np.sin(theta)],
				   [  0                        ,     1                        ,     0            ],
				   [  np.sin(theta)            ,     0                        ,     np.cos(theta)]])
	
	T_Z=np.matrix([[  np.cos(psi)              ,     np.sin(psi)              ,     0            ],
				   [ -np.sin(psi)              ,     np.cos(psi)              ,     0            ],
				   [  0                        ,     0                        ,     1            ]])
	T_BE=T_X*T_Y*T_Z
	
	#Rotation to go from Body Frame to Mars Centered
	T_Cb=np.transpose(T_EC)*np.transpose(T_BE)
	return T_Cb


def T_Eb(theta,phi,psi):
    
	   
	#Rotations to go from Vehicle Carried Normal frame to Body frame
	T_X=np.matrix([[  1 					   ,     0                        ,     0            ],
				   [  0                        ,     np.cos(phi)	          ,     np.sin(phi)	 ],
				   [  0                        ,    -np.sin(phi)              ,     np.cos(phi)  ]])
	
	T_Y=np.matrix([[  np.cos(theta)            ,     0                        ,    -np.sin(theta)],
				   [  0                        ,     1                        ,     0            ],
				   [  np.sin(theta)            ,     0                        ,     np.cos(theta)]])
	
	T_Z=np.matrix([[  np.cos(psi)              ,     np.sin(psi)              ,     0            ],
				   [ -np.sin(psi)              ,     np.cos(psi)              ,     0            ],
				   [  0                        ,     0                        ,     1            ]])
	
	T_BE=T_X*T_Y*T_Z
	T_Eb=np.transpose(T_BE)
	return T_Eb


def T_CE(delta,tau):
	"""Transformation Matrix to go from Vehicle-carried normal Mars frame (FE) to inertial (FI) reference frame
	   delta=Latitude, tau=Longitude
    """
	#Rotation needed to go from Mars Centered to Vehicle Carried Normal Frame
	T_EC=np.matrix([[-np.sin(delta)*np.cos(tau),	-np.sin(delta)*np.sin(tau),		np.cos(delta)],
				    [-np.sin(tau)              ,	 np.cos(tau)              ,     0            ],
					[-np.cos(delta)*np.cos(tau),    -np.cos(delta)*np.sin(tau),    -np.sin(delta)]])
	
	T_CE=np.transpose(T_EC)
	return T_CE

def T_EC(delta,tau):
	"""Transformation Matrix to go from inertial (FI) reference frame to Vehicle-carried normal Mars frame (FE)
	   delta=Latitude, tau=Longitude
    """
	#Rotation needed to go from Mars Centered to Vehicle Carried Normal Frame
	T_EC=np.matrix([[-np.sin(delta)*np.cos(tau),	-np.sin(delta)*np.sin(tau),		np.cos(delta)],
				    [-np.sin(tau)              ,	 np.cos(tau)              ,     0            ],
					[-np.cos(delta)*np.cos(tau),    -np.cos(delta)*np.sin(tau),    -np.sin(delta)]])
	
	return T_EC

def accel(M,Fext,Omega_CI,V,R):
	"""Computes acceleration of translational motion. All inputs except Fext are lists. Omega, V, R are lists of the following form [a,b,c] for given i.
	"""
	Omega_CI=np.array(Omega_CI).reshape(3,1)
	V=np.array(V).reshape(3,1)
	R=np.array(R).reshape(3,1)
	
	acc = 1/M * Fext - 2 * np.cross(Omega_CI,V,axis=0).reshape(3,1) - np.cross(Omega_CI,np.cross(Omega_CI,R,axis=0),axis=0).reshape(3,1)
	return acc

def accel_Omega(I,Mext,dIdt,Omega_bC):
	""" Computes the angular acceleration. I is the inertia tensor (3x3 Matrix), Mext is an array, dIdt is the derivative of the inertia tensor (3x3 matrix), Omega is a list of length 3 containing angular rates
	"""
	Omega_bC=np.array(Omega_bC).reshape(3,1)
	acc_Omega_bC=I.I*(Mext-dIdt*Omega_bC-np.cross(Omega_bC,I*Omega_bC,axis=0))
	return acc_Omega_bC


def get_Vn_dot(M, Fx, Ve, Vn, Vd, R, delta, Omegat):
    Vn_dot=Fx/M - 2 * Omegat * Ve * np.sin(delta) - Omegat * Omegat * R * np.sin(delta) * np.cos(delta) - (Ve * Ve * np.tan(delta) - Vn * Vd)/R
    return Vn_dot

def get_Ve_dot(M, Fy, Ve, Vn, Vd, R, delta, Omegat):
    Ve_dot=Fy/M + 2 * Omegat * ( Vd * np.cos(delta) + Vn * np.sin(delta) )+ Ve/R * (Vn * np.tan(delta) + Vd)
    return Ve_dot

def get_Vd_dot(M, Fz, Ve, Vn, R, delta, Omegat):
    Vd_dot=Fz/M - 2 * Omegat * Ve * np.cos(delta) - Omegat * Omegat * R * np.cos(delta) * np.cos(delta) - (Ve * Ve + Vn * Vn)/R
    return Vd_dot

def get_p_tilde(p, theta, psi, delta_dot, delta, tau_dot, Omegat):
    p_tilde = p +  np.cos(theta) * np.sin(psi) * delta_dot - ( np.cos(delta) * np.cos(psi) * np.cos(theta) + np.sin(delta) * np.sin(theta) ) * (tau_dot + Omegat)
    return p_tilde

def get_q_tilde(q, psi, theta, phi, delta_dot, delta, tau_dot, Omegat):
	q_tilde = q + ( np.sin(psi) * np.sin(theta) * np.sin(phi) + np.cos(psi) * np.cos(phi)) * delta_dot - ( np.cos(delta) * ( np.cos(psi) * np.sin(theta) * np.sin(phi) - np.sin(psi) * np.cos(phi) ) - np.sin(delta) * np.cos(theta) * np.sin(phi) ) * (tau_dot * Omegat)
	return q_tilde

def get_r_tilde(r, psi, phi, theta, delta_dot, delta, tau_dot, Omegat):
	r_tilde = r + ( np.sin(psi) * np.cos(phi) * np.sin(theta) - np.cos(psi) * np.sin(phi) ) * delta_dot - ( np.cos(delta) * ( np.sin(psi) * np.sin(phi) + np.cos(psi) * np.sin(theta) * np.cos(phi) ) - np.sin(delta) * np.cos(theta) * np.cos(phi) ) * (tau_dot * Omegat)
	return r_tilde

def get_phi_dot(p_tilde, q_tilde, r_tilde, phi, theta):
	phi_dot = p_tilde + np.sin(phi) * np.tan(theta) * q_tilde + np.cos(phi) * np.tan(theta) * r_tilde
	return phi_dot

def get_theta_dot(q_tilde, r_tilde, phi):
	theta_dot = np.cos(phi) * q_tilde - np.sin(phi) * r_tilde
	return theta_dot

def get_psi_dot(q_tilde, r_tilde, phi, theta):
	psi_dot = np.sin(phi) / np.cos(theta) * q_tilde + np.cos(phi) / np.cos(theta) * r_tilde
	return psi_dot
#========================================================================================================================================================================================================================================================
#   Input Parameters
#========================================================================================================================================================================================================================================================
################################
# related to Mars & Base       #
################################
Req=3396.2*10**3               #[m] equatorial radius
Rmars=3389.5*10**3             #[m] volumetric mean radius
mu=0.042828*10**6*(10**3)**3   #[m^3/s^2] gravitational parameter  
J2=1960.45*10**(-6)            #[-] scaling factor of the gravitypotential to account for the flattening of Mars
omega_mars=0.7088218*10**(-4)  #[rad/s] Martian angular velocity
M_atm=43.34                    #[g/mol] Mean molecular weight
Rgas=8314.4621/M_atm           #[J/(kg*K)] Martian air gas constant
gamma_gas=1.37                 #[-] heat capacity of Martian air
delta0=41*np.pi/180             #[rad] latitude of the Martian base 
tau0=(23.5)*np.pi/180             #[rad] east longitude of the Martian base
h0=-3*10**3					   #[m] altitude wrt R of Mars base
h_phasing=609.74*10**3         #[m] altitude of the phasing orbit
V_phasing=3272.466             #[m/s]
################################
# related to the vehicle 	   #
################################
tb=489.5180709799191		   #[s] burn time of the engines 489.5180709799191
initial_tilt=(0.51)*np.pi/180      #[rad] initial tilt off the vertical axis 0.51*np.pi/180
#Propulsive Parameters
Isp= 383.250565907662          #[s] LOX-LCH4 ##CAN BE UPDATED##
ceff=get_ceff(Isp)             #[m/s]
n=9                            #[-] number of engines             ##CAN BE UPDATED##
De=1.35049466031671            #[m] Diameter of an engine         ##CAN BE UPDATED##
Ae=np.pi/4*De*De               #[m] Exit area of engine
pe=6077.910186177842           #[Pa] exit pressure of exhaust gasses ##CAN BE UPDATED##
#Structural Parameters
d=6                            #[m] diameter of the vehicle assuming cylinder ##CAN BE UPDATED##
lcap=5.16	                   #[m] length of the capsule          ##CAN BE UPDATED##
Acap=22.144                    #[m^2] lateral cross-sectional area of capsule (Assuming parabolic form) ##CAN BE UPDATED##
LDratio=3.2                    #[-] Length-over-Diameter ratio     ##CAN BE UPDATED##
ltot=LDratio*d                 #[m] overall length of the vehicle 
leng=3  	                   #[m] lenght of engine               ##CAN BE UPDATED##
lcyl=ltot-lcap-leng            #[m] length of everything below the capsule
Acyl=lcyl*d                    #[m^2] lateral cross-sectional area of cylyinder below capsule
xcg0=5.309139075 			   #[m] c.g. position on ground from base of the cylinder ##CAN BE UPDATED##
t_xcg0=0                       #[s] time at which xcg0 happens
xcgm=5.379751314               #[m] c.g. position at mdot max point ##CAN BE UPDATED##
# print this to get the time at which mdot happens others["time"][list(others["mdot"]).index(np.max(others["mdot"]))]
t_xcgm=131.78999999991152      #[s] time at which xcg of max fuel flow happens
xcge=5.400472234               #[m] c.g. position at end of ascent from base of the cylinder ##CAN BE UPDATED##
t_xcge=tb
Iy0=8945330.201                #[kg*m^2] MMOI x on ground
Ix0=844679.9158				   #[kg*m^2] MMOI y on ground
Iz0=8945330.201                #[kg*m^2] MMOI z on ground
Iym=6794407.105	               #[kg*m^2] MMOI x point of mdot max
Ixm=802779.394			       #[kg*m^2] MMOI y point of mdot max
Izm=6794407.105                #[kg*m^2] MMOI z point of mdot max
Iye=4095074.195                #[kg*m^2] MMOI x at burn out
Ixe=328306.1762				   #[kg*m^2] MMOI y at burn out
Ize=4095074.195                #[kg*m^2] MMOI z at burn out
#Aerodynamic Parameters
S=np.pi/4*d**2                 #[m^2] reference surface area used here is the area of the circle of the cross-section  
#========================================================================================================================================================================================================================================================
#   Values needed by others to update some Input Parameters
#========================================================================================================================================================================================================================================================
# To update MOIs and xcg:
if updateMOI:
	others = ascent_sim(tb=148.7274456555216,initial_tilt=3.2,i_base=42.5,h0=-3*10**3,d=7.67,M_initial=198948.0,Mp_class2=140235.07208472778,Isp=350,n=6,De=2,pe=5066.25)[4]
	Mp1,Mp2,MpP1,MpP2=inputDimitris(others["Mprop"],162719,others["mdot"])
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("Point between zero and mdot max: ",MpP1)
	print("Point of mdot max.: ",Mp1)
	print("Point between mdot max and mdot final: ",MpP2)
	print("Point of mdot final: ",Mp2)
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	
#========================================================================================================================================================================================================================================================
#   Initial Conditions  
#========================================================================================================================================================================================================================================================
#load the pitch rate due to the gravity gotten from the simplified ascent simulation
gamma_dot=np.genfromtxt("gravity_pitchrate.txt")
t=[0]					       #[s] time of flight 
dt=0.1					       #[s] step size
#Mass
M=[202413.4011]
r0=Rmars+h0+100                     
R=[r0]
Vn=[0]
Ve=[0]
Vd=[110/3.6]
Vnorm=[np.sqrt(Vn[0]*Vn[0]+Vd[0]*Vd[0]+Ve[0]*Ve[0])]
delta=[delta0]
tau=[tau0]
psi=[np.pi/2]
theta=[-(np.pi/2-initial_tilt)]
phi=[0]
p_dot=[0]
q_dot=[0]
r_dot=[0]
p=[0]
q=[0]
r=[0]
#am=aerodynamic moment to be coutneracted by control to be stable
am=[]
#prescribed max TWratio
TW0=1.5
TWe=4
mdot_list=[]
ac=[]
Fthrust=[]
deltaYaw=[]
xcp=[get_xcp(lcap,lcyl,ltot,xcg0,Acap,Acyl,delta_cg(t[-1],xcg0,t_xcg0,xcgm,t_xcgm,xcge,t_xcge,tb))]          #[m] x-location of center of pressure from the c.g. (i.e. the origin)
#========================================================================================================================================================================================================================================================
#   Simulation
#========================================================================================================================================================================================================================================================
i=0
while R[-1]<Rmars+h_phasing:
#for i in range(2):
	i+=1
	########################################
	# Compute atmospheric properties at i  #
	########################################
	T=get_T(R[-1]-Rmars)
	pr=get_p(R[-1]-Rmars)
	rho=get_rho(pr,T,Rgas)
	a=get_a(gamma_gas,Rgas,T)
	Mach=get_Mach(Vnorm[-1],a)
	Cl,Cd=SS_aerodynamics_coefficients(Mach,0)
	qdyn=0.5*rho*Vnorm[-1]*Vnorm[-1]						#dynamic pressure
	g=get_g(mu,Req,Rmars,R[-1]-Rmars,delta,J2)
	#magnitude of the thrust force
	if t[-1]>=tb:
		Ftmag=0
		mdot=0
	else:
		TWratio=TWlinear(t[-1],0,tb,TW0,TWe)
		Ftmag=(TWratio*M[-1]*g[-1])
		mdot=(TWratio*M[-1]*g[-1]-Ae*(pe-pr))/ceff
	Fthrust.append(Ftmag)
	mdot_list.append(mdot)
    ########################################
    # Compute external Forces              #
    ########################################
    #external forces expressed in the vehicle normal carried reference frame
	Fext = T_Eb(theta[-1],phi[-1],psi[-1]) * np.array([[ -Ftmag - M[-1] * g[-1] * np.sin(theta[-1]) - qdyn * S * Cd[-1]           ],
	                                                   [ M[-1] * g[-1] * np.cos(theta[-1]) * np.sin(phi[-1])                     ],
	                                                   [ M[-1] * g[-1] * np.cos(theta[-1]) * np.cos(phi[-1]) - qdyn * S * Cl[-1] ]])
	Fx=float(Fext[0])
	Fy=float(Fext[1])
	Fz=float(Fext[2])
    ########################################
    # Compute external Moments             #
    ########################################
    #in the body frame!
	Mcm =  0 * qdyn * S * Cl[-1] * xcp[-1]
	am.append(qdyn * S * Cl[-1] * xcp[-1])
    ########################################
	# Compute Accelerations                #
	########################################	
	Vn_dot=get_Vn_dot(M[-1], Fx, Ve[-1], Vn[-1], Vd[-1], R[-1], delta[-1], Omegat=omega_mars)
	Ve_dot=get_Ve_dot(M[-1], Fy, Ve[-1], Vn[-1], Vd[-1], R[-1], delta[-1], Omegat=omega_mars)
	Vd_dot=get_Vd_dot(M[-1], Fz, Ve[-1], Vn[-1], R[-1], delta[-1], Omegat=omega_mars)
	#break
	ac.append(np.sqrt(Vn_dot*Vn_dot + Ve_dot*Ve_dot + Vd_dot*Vd_dot))
    ########################################
	# Compute Velocities                   #
	########################################
	Vn_new = Vn[-1] + dt * Vn_dot
	Ve_new = Ve[-1] + dt * Ve_dot
	Vd_new = Vd[-1] + dt * Vd_dot
	Vnorm.append(np.sqrt(Vn_new*Vn_new + Ve_new*Ve_new + Vd_new*Vd_new)) 
    ########################################
	# Compute Positions                    #
	########################################
	R_dot = - Vd[-1]
	delta_dot = Vn[-1] / R[-1]
	tau_dot = Ve[-1] / ( R[-1] * np.cos(delta[-1]))
	
	Rnew = R[-1] + dt * R_dot
	deltanew = delta[-1] + dt * delta_dot
	taunew = tau[-1] + dt * tau_dot
	
	print()
	print("Altitude is: ", Rnew-Rmars)
	print("time is: ",t[-1])
	
    ########################################
	# Compute MOIs                         #
	########################################
	if t[-1]<=tb:
	    Ix=get_Ix(t[-1],Ix0,t_xcg0,Ixm,t_xcgm,Ixe,t_xcge,tb)[0] 
	    Iy=get_Iy(t[-1],Iy0,t_xcg0,Iym,t_xcgm,Iye,t_xcge,tb)[0]
	    Iz=get_Iz(t[-1],Iz0,t_xcg0,Izm,t_xcgm,Ize,t_xcge,tb)[0]
	#print(Iz)
	I=np.matrix([[Ix, 0 , 0 ],
			     [0 , Iy, 0 ],
				 [0 , 0 , Iz]])
    ########################################
	# Compute derivative of MOIs           #
	########################################
	if t[-1]<=tb:
	    Ixdot=get_Ix(t[-1],Ix0,t_xcg0,Ixm,t_xcgm,Ixe,t_xcge,tb)[1] 
	    Iydot=get_Iy(t[-1],Iy0,t_xcg0,Iym,t_xcgm,Iye,t_xcge,tb)[1]
	    Izdot=get_Iz(t[-1],Iz0,t_xcg0,Izm,t_xcgm,Ize,t_xcge,tb)[1]
	#print(Iz)
	dIdt=np.matrix([[Ixdot, 0    , 0    ],
			        [0    , Iydot, 0    ],
				    [0    , 0    , Izdot]])
    ########################################
	# Compute angular accelerations & rates#
	########################################
	Omega_bI=np.array([ [ p[-1] ], [ q[-1] ], [ r[-1] ] ])		
	Omega_bI_dot = I.I * ( - dIdt * Omega_bI - np.cross( Omega_bI, I * Omega_bI ,axis=0) )
	
	# ?? add control over pitch and yaw ??
	
	p_dotnew=float(Omega_bI_dot[0])
	q_dotnew=float(Omega_bI_dot[1])
	r_dotnew=float(Omega_bI_dot[2])
	p_dot.append(p_dotnew)
	q_dot.append(q_dotnew)
	r_dot.append(r_dotnew)
	
	pnew = p[-1] + dt * p_dotnew
	#qnew = q[-1] + dt * q_dotnew
	qnew=gamma_dot[i]
	rnew = r[-1] + dt * r_dotnew
	
    ########################################
	# Compute Euler angles                 #
	########################################	
	p_tilde=get_p_tilde(p[-1], theta[-1], psi[-1], delta_dot, delta[-1], tau_dot, Omegat=omega_mars)
	q_tilde=get_q_tilde(q[-1], psi[-1], theta[-1], phi[-1], delta_dot, delta[-1], tau_dot, Omegat=omega_mars)
	r_tilde=get_r_tilde(r[-1], psi[-1], phi[-1], theta[-1], delta_dot, delta[-1], tau_dot, Omegat=omega_mars)
	
	phi_dot=get_phi_dot(p_tilde, q_tilde, r_tilde, phi[-1], theta[-1])
	theta_dot=get_theta_dot(q_tilde, r_tilde, phi[-1])
	psi_dot=get_psi_dot(q_tilde, r_tilde, phi[-1], theta[-1])
	
	phinew = phi[-1] + dt * phi_dot
	thetanew = theta[-1] + dt * theta_dot
	psinew = psi[-1] + dt * psi_dot
	
	
	
	#append velocities
	Vn.append(Vn_new)
	Ve.append(Ve_new)
	Vd.append(Vd_new)
	#append delta
	delta.append(deltanew)
	#append tau
	tau.append(taunew)
	#append R
	R.append(Rnew)
	#append angular rates
	p.append(pnew)
	q.append(qnew)
	r.append(rnew)
    #append Euler angles
	phi.append(phinew)
	theta.append(thetanew)
	psi.append(psinew)
    ########################################
    # Compute new Mass, c.p. & angles      #
    ########################################
	M.append( M[-1] - mdot * dt)
    #update xcp
	xcp.append(get_xcp(lcap,lcyl,ltot,xcg0,Acap,Acyl,delta_cg(t[-1],xcg0,t_xcg0,xcgm,t_xcgm,xcge,t_xcge,tb)))
    ########################################
    # Update the time                      #
    ########################################
	t.append(t[-1] + dt)
	
#========================================================================================================================================================================================================================================================
#   Reformatting
#========================================================================================================================================================================================================================================================
#total propellant mass needed
Mp=np.sum(np.array(mdot_list)*dt)
#delta V needed
DeltaV = ceff * np.log( M[0] / ( M[0] - Mp ) )


M=np.array(M)
r0=Rmars+h0+100                     
R=np.array(R)
V=np.array(Vnorm)
delta=np.array(delta)
tau=np.array(tau)
psi=np.array(psi)
theta=np.array(theta)
phi=np.array(psi)
p_dot=np.array(p_dot)
q_dot=np.array(q_dot)
r_dot=np.array(r_dot)
p=np.array(p)
q=np.array(q)
r=np.array(r)
Ft=np.array(Fthrust)
#am=aerodynamic moment to be coutneracted by control to be stable
am=np.array(am)
#prescribed max TWratio
TW0=1.5
TWe=4
mdot=np.array(mdot_list)
a_array=np.array(ac)			 
t=np.array(t)	
				 
#compute cartesian coordinates of the trajectory
X = R * np.cos(delta) * np.cos(tau)
Y = R * np.cos(delta) * np.sin(tau)
Z = R * np.sin(delta)

#Save the pitch acceleration 
np.savetxt("pitch_angular_acceleration.txt", list(q_dot), delimiter=",")
#========================================================================================================================================================================================================================================================
#   Plotting
#========================================================================================================================================================================================================================================================

if plotting:
	
	########################################
	# 2D plots                             #
	########################################
	#t vs R
    plt.figure()
    plt.plot(t,(R-Rmars)/10**3,color="navy")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Altitude")
    plt.xlabel("Time [s]")
    plt.ylabel("Altitude [km]")
    plt.show()

    #t vs Ft
    ac_max=np.ones((len(t),1))*4
    plt.figure()
    plt.plot(t[1::],a_array*3.71/9.81,color="hotpink")
    plt.plot(t,ac_max,linestyle=":",color="firebrick")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Acceleration")
    plt.xlabel("Time [s]")
    plt.ylabel("Acceleration [g_earth's]")
    plt.show()


    #t vs V
    Vreq=np.ones((len(Vnorm),1))*V_phasing
    plt.figure()
    plt.plot(t,V/10**3,color="blue")
    plt.plot(t,Vreq/10**3,linestyle=":",color="firebrick")    
    plt.grid(color="gainsboro")
    plt.title("Time vs Velocity")
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    plt.show()
		
    #t vs Thrust
    plt.figure()
    plt.plot(t[1::],Ft/10**3,color="maroon") 
    plt.grid(color="gainsboro")
    plt.title("Time vs Thrust")
    plt.xlabel("Time [s]")
    plt.ylabel("Thrust [kN]")
    plt.show()

    #t vs Euler angles
    plt.figure()
    plt.plot(t,phi,color="red", label="Phi (Roll Angle)") 
    plt.plot(t,theta,color="blue", label="Theta (Pitch Angle)")
    plt.plot(t,psi,color="green",  label="Psi (Yaw Angle)")
    plt.grid(color="gainsboro")
    plt.title("Time vs Pitch Angle")
    plt.xlabel("Time [s]")
    plt.ylabel("Angles [deg]")
    plt.legend()
    plt.show()


	
	#Plot for visualizing how MMOI and cg shift along flight. Also, the parabolic fit is added.
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('MMOI [kg*m^2]', color="tab:red")
    ax1.plot([0,t_xcgm,tb], [Ix0,Ixm,Ixe], color="orange", label="I_x")
    ax1.plot(np.linspace(0,tb,1000), get_Ix(np.linspace(0,tb,1000),Ix0,t_xcg0,Ixm,t_xcgm,Ixe,t_xcge,tb)[0], linestyle=":" , color="moccasin", label="I_x: parabolic fit")
    ax1.plot([0,t_xcgm,tb], [Iy0,Iym,Iye], color="tab:red", label="I_y")
    ax1.plot(np.linspace(0,tb,1000), get_Iy(np.linspace(0,tb,1000),Iy0,t_xcg0,Iym,t_xcgm,Iye,t_xcge,tb)[0], linestyle=":" , color="indianred", label="I_y: parabolic fit")
    ax1.plot([0,t_xcgm,tb], [Iz0,Izm,Ize], color="tab:red", label="I_z")
    ax1.plot(np.linspace(0,tb,1000), get_Iz(np.linspace(0,tb,1000),Iz0,t_xcg0,Izm,t_xcgm,Ize,t_xcge,tb)[0], linestyle=":" , color="indianred", label="I_z: parabolic fit")
    ax1.tick_params(axis='y', labelcolor="tab:red")
    plt.legend(loc="upper left")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis


    ax2.set_ylabel('X_c.g. [m]', color="tab:blue")
    ax2.plot([0,t_xcgm,tb], [xcg0,xcgm,xcge], color="tab:blue", label="X_c.g.")
    ax2.plot(np.linspace(0,tb,1000), xcg(np.linspace(0,tb,1000),xcg0,t_xcg0,xcgm,t_xcgm,xcge,t_xcge,tb), linestyle=":" , color="tab:blue", label="X_c.g.: parabolic fit")
    ax2.tick_params(axis='y', labelcolor="tab:blue")
    plt.grid(color="gainsboro")
    plt.title("MMOI and X_c.g. variation in time")
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend(loc="upper right")
    plt.show()
 

   
    ########################################
    # 3D plots                             #
    ########################################
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Make data
    r = Rmars
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    x1 = r*sin(phi)*cos(theta)
    y1 = r*sin(phi)*sin(theta)
    z1 = r*cos(phi)


    # Plot the surface
    plt.gca().patch.set_facecolor('white')
    ax.w_xaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
    ax.w_yaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
    ax.w_zaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
    ax.plot_surface(x1, y1, z1, color='coral')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')    
    ax.scatter(X[0],Y[0],Z[0],color="lime")


    #Plot of trajectory X,Y,Z
    ax.plot(X,Y,Z,color="dodgerblue")
    plt.show()

    #Plot trajectory on its own
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(X,Y,Z,color="dodgerblue")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z') 
    plt.show()


print("[---", (time.time() - start_time) ,"seconds ---]" )