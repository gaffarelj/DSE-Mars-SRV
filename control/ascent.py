import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
from reentry_footprint import flight, time, dt, mean_radius, mars, rho_0
import disturbances as dist

import numpy as np
from matplotlib import pyplot as plt

#Initial slew values
#Assume spin acceleration/deceleration of 5%, coast time of 90%
angle = 30 * np.pi / 180
slew_duration = 100 #s
slew_acc_duration = 0.05 * slew_duration

def slew_ascent(angle,slew_duration):
    spin_rate =  slew_angle_tot / slew_duration
    spin_acc  =  spin_rate      / slew_acc_duration
    spin_dec  = -spin_acc
angle = 30
slew_duration = 100 #s
slew_acc_duration = 0.05 * slew_duration
