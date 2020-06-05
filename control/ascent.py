import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act

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
    Mp = act.RCSpropellant(RCS_impulse,slew_duration,Isp)
    return RCS_thrust, RCS_impulse, Mp

RCS_thrust, RCS_impulse, Mp = slew_ascent(angle,slew_duration,Iy,cg_orbit)
RCS_thrust_max = np.max(RCS_thrust)
print('thrust per RCS engine:',RCS_thrust_max)
print('Impulse per RCS engine:', RCS_impulse)
print('Total propellant needed:', Mp)
