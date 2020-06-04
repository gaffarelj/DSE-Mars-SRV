import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT


def get_h(planet, V):
	return planet.mu / V ** 2 - planet.r

def get_V(planet, h):
	return np.sqrt(planet.mu / (planet.r + h))

mars = AT.Planet()
# 175 km --> 108.1 hrs (4.4 days) --> 4 days
h0 = 175*1e3
he = 125*1e3
V0 = get_V(mars, h0)
#print(V0, V0-1427.87)
S = 11.9
m = 12500 - 1500
gamma = 0

mars = AT.Planet()
state = np.zeros(6)
state[0] = V0						# velocity
state[1] = np.radians(gamma)		# flight path angle
state[2] = np.radians(gamma)		# heading angle
state[3] = mars.r + h0				# radius
state[4] = np.radians(-27.088)      # longitude 25.5
state[5] = np.radians(4.51)         # latitude 42.5
roll_angle = np.radians(0)
aoa = np.radians(0)

cd = 1.2
cl = 0.23
S = 11.9

motion = AT.Motion(state, roll_angle, aoa, S, m, cl, cd, mars, end_t=7 * 24 * 3600)
flight, time = motion.forward_euler(0.1)

t, a, v, q, h, angle, m = time, np.array(motion.a_s), flight[:,0], np.array(motion.q_s), flight[:,3] - mars.r, flight[:,1], motion.mass

AT.plot_dual([tt / 60 / 60 for tt in t], v, [ht / 1000 for ht in h], "Time [hr]", "Velocity [m/s]", "Altitude [km]")
AT.plot_dual([tt / 60 / 60 for tt in t], q, [ht / 1000 for ht in h], "Time [hr]", "Dyn. press [Pa]", "Altitude [km]")

try:
	reenter_i = np.where(h <= he)[0][0]
except IndexError:
	reenter_i = -1
reenter_t = t[reenter_i]

print(f"Starting from {h0//1000} km ({int(V0)} m/s), below h={he//1000} km after {round(reenter_t/3600, 1)} hrs")