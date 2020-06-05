import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act
import reentry_control as rty
import rendez_vous as rvs

import numpy as np
from matplotlib import pyplot as plt

#Initial slew values
#Assume spin acceleration/deceleration of 5%, coast time of 90%
angle = 20 * np.pi / 180
slew_duration = 20 #s

#Vehicle constants
length = act.length_body
width  = act.body_radius * 2
Ix = act.Ix
Iy = act.Iy
Iz = act.Iz
cg_orbit = act.z_cg_orbit
error_angle = 2 * np.pi / 180

#propellant properties
Isp = act.Isp


def slew_ascent(slew_angle,slew_duration,I,cg):
    slew_acc_duration = 0.05 * slew_duration
    slew_dec_duration = slew_acc_duration
    spin_rate =  slew_angle / slew_duration
    spin_acc  =  spin_rate  / slew_acc_duration
    spin_dec  = -spin_acc

    RCS_torque = (I * spin_acc)
    RCS_thrust = act.RCS_torque_to_thrust(RCS_torque,"y",cg,'normal')
    RCS_impulse =  RCS_thrust * (slew_acc_duration + slew_dec_duration)
    Mp = RCS_impulse / (act.Isp * 9.80665)
    return RCS_torque, RCS_thrust, RCS_impulse, Mp

RCS_torque, RCS_thrust, RCS_impulse, Mp = slew_ascent(angle,slew_duration,Iy,cg_orbit)
RCS_thrust_max = np.max(RCS_thrust)

T_error_z, T_error_y, T_error_x = act.thrust_error(RCS_thrust,cg_orbit,error_angle)
RCS_error_x  = act.RCS_torque_to_thrust(T_error_x,'y',cg_orbit,'error_bottom')
RCS_error_y  = act.RCS_torque_to_thrust(T_error_y,'y',cg_orbit,'error_bottom')
RCS_error_z  = act.RCS_torque_to_thrust(T_error_z,'y',cg_orbit,'error_bottom')
RCS_error    = max([RCS_error_x,RCS_error_y,RCS_error_z])
mp_error     = act.RCSpropellant(RCS_error,slew_duration,Isp)

RCS_failure  = act.RCS_torque_to_thrust(RCS_torque,'y',cg_orbit,'failure') - RCS_thrust
# mp_failure   = Mp * (RCS_failure/RCS_thrust)

print('thrust per engine:',RCS_thrust_max)
print('Impulse per RCS engine:', RCS_impulse)
print('Total propellant needed:', Mp)

print('Redundancy thrust per engine:', RCS_error + RCS_failure)
print('Redundancy propellant:', mp_error)
print('Total redundant propellant needed:', Mp + mp_error)


Mpropellant_total = rvs.mp_total + Mp + mp_error + 95.51959871385976 + rty.mp_roll

print(Mpropellant_total)
