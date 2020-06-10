import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act
import reentry_control as rty
import rendez_vous as rvs

import numpy as np
from matplotlib import pyplot as plt

margin = 2.



#=====================================================================================================================================================================================================================
#Node properties
#=====================================================================================================================================================================================================================
mu     = 0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
R      = 3389.5*10**3                  #[m] volumetric mean radius
h_node = 500*10**3
h_phasing=609.74*10**3
g_earth= 9.81             #[m]

period =  2*np.pi* np.sqrt((R+h_node+h_phasing) ** 3 / mu)
omega  = (2 * np.pi)/period
V      = omega * (R+h_node+h_phasing)

#=====================================================================================================================================================================================================================
#Initial slew and Delta V values
#=====================================================================================================================================================================================================================
#Assume spin acceleration/deceleration of 5%, coast time of 90%
angle = (180. - 26.7) * np.pi / 180
slew_duration = 0.5 * period #s

# slew_duration2 = 121.
# deltaV = 0.014 * V

#=====================================================================================================================================================================================================================
#Vehicle constants
#=====================================================================================================================================================================================================================
length = act.length_body
width  = act.body_radius * 2
Ix = 4098123.837
Iy = 315781.2849
Iz = 4098123.837
cg_orbit = act.z_cg_orbit
m  = 53056.28092

#propellant properties
Isp_mono = 140
Isp      = 317

#=====================================================================================================================================================================================================================
#Function to compute thrust from slew maneuver
#=====================================================================================================================================================================================================================
def slew_ascent(slew_angle,slew_duration,I,cg,Isp):
    slew_acc_duration = 0.5 * slew_duration
    slew_dec_duration = slew_acc_duration
    spin_rate =  slew_angle / slew_duration
    spin_acc  =  spin_rate  / slew_acc_duration
    spin_dec  = -spin_acc

    RCS_torque = (I * spin_acc)
    RCS_thrust = margin * act.RCS_torque_to_thrust(RCS_torque,"y",cg,'normal')
    RCS_impulse =  RCS_thrust * (slew_acc_duration + slew_dec_duration)
    Mp = RCS_impulse / (Isp * 9.80665)
    return RCS_torque, RCS_thrust, RCS_impulse, Mp

#=====================================================================================================================================================================================================================
# Function to compute thrust from Delta V
#=====================================================================================================================================================================================================================
def vac_thrust(DeltaV,Isp,Mbegin,tb,De=0,pe=0):
    """computes vacuum thrust from DeltaV.
        Mbegin=mass at the start of the maneuver, tb=burn time, De=exit diameter of engine/thruster, pe=exhaust exit pressure of engine/thruster
    """
    Ae=np.pi/4*De*De
    thrust=(np.exp(DeltaV/(Isp*9.80665))*Mbegin-Mbegin)/(tb*np.exp(DeltaV/(Isp*9.80665)))*Isp*9.80665+Ae*(pe)
    Mprop=(np.exp(DeltaV/(Isp*9.80665))*Mbegin-Mbegin)/(np.exp(DeltaV/(Isp*9.80665)))

    return thrust, Mprop


#========================================================================
# Slew maneuver
#========================================================================
RCS_torque, RCS_thrust, RCS_impulse, Mp = slew_ascent(angle,slew_duration,Iy,cg_orbit,Isp_mono)


#=====================================================================================================================================================================================================================
#Errors
#=====================================================================================================================================================================================================================
#========================================================================
# Thruster misalignment
#========================================================================
error_angle = 2 * np.pi / 180

T_error_y, T_error_x, T_error_z = act.thrust_error(RCS_thrust,cg_orbit,error_angle)
#========================================================================
# Disturbances
#========================================================================
theta = 2. * np.pi / 180

Tgx  = dist.gravitygradient_disturbance(Iy,Iz,omega,theta)
Tgy  = dist.gravitygradient_disturbance(Ix,Iz,omega,theta)
Tsp = dist.solarpressure_disturbance(theta,cg_orbit)
Tm  = dist.magnetic_disturbance(R)

#========================================================================
# Total
#========================================================================
T_error_x += Tgx + Tsp + Tm
T_error_y += Tgy
RCS_error_x  = act.RCS_torque_to_thrust(T_error_x,'y',cg_orbit,'error_bottom')
RCS_error_y  = act.RCS_torque_to_thrust(T_error_y,'x',cg_orbit,'error_bottom')
RCS_error_z  = act.RCS_torque_to_thrust(T_error_z,'z',cg_orbit,'error_bottom')

RCS_error    = max([RCS_error_x,RCS_error_y,RCS_error_z])
mp_error     = act.RCSpropellant(RCS_error,slew_duration,Isp)

# RCS_failure  = act.RCS_torque_to_thrust(RCS_torque,'y',cg_orbit,'failure') - RCS_thrust
# mp_failure   = Mp * (RCS_failure/RCS_thrust)
print('thrust per engine:',RCS_torque)
print('thrust per engine:',RCS_thrust)
print('Impulse per RCS engine:', RCS_impulse)
print('Total propellant needed:', Mp + mp_error)
print('REDUNDANCY')
print('Misalignment torque: ', T_error_x-Tgx-Tsp-Tm,T_error_y-Tgy,T_error_z)
print('Disturbance torque: ', Tgx+Tsp+Tm,Tgy,0)
print('Redundancy thrust per engine:', RCS_error)
print('Redundancy propellant:', mp_error)
print('Total redundant propellant needed:', mp_error)


Mpropellant_total = rvs.mp_total + Mp + mp_error + rty.mp
print(Mpropellant_total)
