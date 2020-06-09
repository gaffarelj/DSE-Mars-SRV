import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT
from astrodynamics import BAT


def run_motion(V0, gamma, h0, chutes=[], print_deploy=False, prop_reentry=[], back_orbit=[]):
	mars = AT.Planet()
	state = np.zeros(6)
	state[0] = V0						# velocity
	state[1] = np.radians(gamma)		# flight path angle
	state[2] = np.radians(gamma)		# heading angle
	state[3] = mars.r + h0				# radius 
	state[4] = np.radians(-27.088)      # longitude 25.5
	state[5] = np.radians(4.51)         # latitude 42.5

	roll_angle = np.radians(0)
	aoa = np.radians(-5)

	cd = 1.2
	cl = 0.23
	S = 11.9

	vehicle_mass = 14200-1203

	motion = AT.Motion(state, roll_angle, aoa, S, vehicle_mass, cl, cd, mars, chutes, print_deploy, prop_reentry, end_t=5*24*3600, back_orbit=back_orbit)
	flight, time = motion.forward_euler(0.1)
	return time, np.array(motion.a_s), flight[:,0], np.array(motion.q_s), flight[:,3] - mars.r, flight[:,1], motion.mass

V0 = 3318-261
gamma=0
h0=500e3
t, a, v, q, h, angle, m = run_motion(V0, gamma, h0, back_orbit=[425, 260e3])
AT.plot_dual([tt/3600 for tt in t], v, [ht/1000 for ht in h], "Time [hr]", "Velocity [m/s]", "Altitude [km]")