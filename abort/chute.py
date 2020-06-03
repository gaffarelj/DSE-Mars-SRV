import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT
from astrodynamics import BAT


ascent = BAT.ascent_sim()[4]
plt.plot(ascent["time"], ascent["q"])
plt.show()
plt.plot(ascent["time"], ascent["V"])
plt.show()


mars = AT.Planet()

h, v = 15, 600

state = np.zeros(6)
state[0] = v						# velocity
state[1] = np.radians(85)			# flight path angle
state[2] = np.radians(85)			# heading angle
state[3] = mars.r + h*10**3			# radius 
state[4] = np.radians(-27.088)      # longitude 25.5
state[5] = np.radians(4.51)         # latitude 42.5

roll_angle = np.radians(0)
aoa = np.radians(5)

cd = 1.2
cl = 0.23
S = 11.9
# q = 50Pa, q ~ 250Pa?
t1, t2 = 282, 315
ballute = AT.pc(0.45, 105, 50, deploy_time=t1, n=1)
rogue = AT.pc(0.6, 95, 120, deploy_time=t1, n=3)
parachute = AT.pc(0.8, 600, 250, deploy_time=t2, n=3)
chutes = [rogue, parachute]

vehicle_mass = 12000
motion = AT.Motion(state, roll_angle, aoa, S, vehicle_mass, cl, cd, mars, parachutes=chutes)
flight, time = motion.forward_euler(0.1)

AT.plot_dual(time, [a/9.81 for a in motion.a_s], flight[:,0], 'Time [s]', 'Acceleration [g]', 'Velocity [m/s]')
AT.plot_dual(time, (flight[:,3] - mars.r)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
AT.plot_dual(time, motion.q_s, (flight[:,3] - mars.r)/1000, 'Time [s]', 'Dyn. Pressure [Pa]', 'Altitude [km]')