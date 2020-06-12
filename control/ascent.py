import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import DYNAS as dns
import disturbances as dist
import actuator_properties as act

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
slew_limit = 0.5 * period #s
tburn = 1.
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
Isp      = 140

#=====================================================================================================================================================================================================================
#Function to compute thrust from slew maneuver
#=====================================================================================================================================================================================================================
def slew_ascent(tburn,slew_angle,I,cg,Isp):
    RCS_thrust = 450.
    torque = act.RCS_thrust_to_torque(RCS_thrust,'z',cg)
    spin_acc   = torque / I
    spin_rate  = spin_acc * tburn
    slew_time = slew_angle / spin_rate
    RCS_impulse = 4 * RCS_thrust * tburn * 2
    Mp = RCS_impulse / (Isp * 9.80665)
    return torque, RCS_thrust, RCS_impulse, Mp,slew_time


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
RCS_torque, RCS_thrust, RCS_impulse, Mp,slew_time = slew_ascent(tburn,angle,Iy,cg_orbit,Isp_mono)
print('Rotation limit: ',slew_limit)
print('Rotation duration:' ,slew_time)
#========================================================================
# Aerodynamic pitch control
#========================================================================
dt = 0.1
am = dns.am
tm = dns.t
pitch_mp = 0
pitch_thrust = act.RCS_torque_to_thrust(tm,'z',cg_orbit,'normal')

for thrust in pitch_thrust:
    pitch_mp += 4 * act.RCSpropellant(thrust,dt,Isp)

#=====================================================================================================================================================================================================================
#Errors
#=====================================================================================================================================================================================================================
#========================================================================
# Thruster misalignment
#========================================================================
error_angle = 2 * np.pi / 180

T_error_y, T_error_x, T_error_z = act.thrust_error(RCS_thrust,cg_orbit,error_angle)
T_error_pitch_y, T_error_pitch_x, T_error_pitch_z = act.thrust_error(max(pitch_thrust),cg_orbit,error_angle)
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
mp_error     = 18 * act.RCSpropellant(RCS_error,tburn,Isp)

T_error_pitch_x += Tgx + Tsp + Tm
T_error_pitch_y += Tgy
RCS_error_pitch_x  = act.RCS_torque_to_thrust(T_error_pitch_x,'y',cg_orbit,'error_bottom')
RCS_error_pitch_y  = act.RCS_torque_to_thrust(T_error_pitch_y,'x',cg_orbit,'error_bottom')
RCS_error_pitch_z  = act.RCS_torque_to_thrust(T_error_pitch_z,'z',cg_orbit,'error_bottom')

RCS_pitch_error  = max([RCS_error_pitch_x,RCS_error_pitch_y,RCS_error_pitch_z])
mp_pitch_error     = 18 * act.RCSpropellant(RCS_pitch_error,tm[-1],Isp)
#=====================================================================================================================================================================================================================
#Total
#=====================================================================================================================================================================================================================
mp_total = Mp + mp_error + pitch_mp + mp_pitch_error
print('Torque:',RCS_torque)
print('Rotation thrust in total  :',4*RCS_thrust)
print('Rotation thrust per engine:',RCS_thrust)
print('Max. pitch control thrust in total:',max(pitch_thrust)*4)
print('Max. pitch control thrust per engine:',max(pitch_thrust))
print('Total propellant needed:', Mp + mp_error+pitch_mp)
print('REDUNDANCY')
print('Misalignment torque rotation     : ', T_error_x-Tgx-Tsp-Tm,T_error_y-Tgy,T_error_z)
print('Misalignment torque pitch control: ', T_error_pitch_x-Tgx-Tsp-Tm,T_error_pitch_y-Tgy,T_error_pitch_z)
print('Disturbance torque: ', Tgx+Tsp+Tm,Tgy,0)
print('Redundancy thrust per engine (rotation):', RCS_error)
print('Redundancy thrust per engine (pitch control):', RCS_pitch_error)

print('Redundancy propellant (rotation):', mp_error)
print('Redundancy propellant (pitch control):', mp_pitch_error)
