import numpy as np
import sys
sys.path.append(".")
from astrodynamics import astro_tools as AT




mars = AT.Planet()

state = np.zeros(6)
state[0] = 400						# velocity
state[1] = np.radians(85)			# flight path angle
state[2] = np.radians(85)			# heading angle
state[3] = mars.r + 10000			# radius 
state[4] = np.radians(-27.088)      # longitude 25.5
state[5] = np.radians(4.51)         # latitude 42.5

roll_angle = np.radians(0)
aoa = np.radians(5)

cd = 1.2
cl = 0.23
S = 11.9

ballute = AT.pc(0.4, 63.6, 25)
parachute = AT.pc(0.8, 600, 50, deploy_time=175, n=3)
chutes = [parachute]

vehicle_mass = 11500
motion = AT.Motion(state, roll_angle, aoa, S, vehicle_mass, cl, cd, mars, parachutes=chutes)
flight, time = motion.forward_euler(0.1)

AT.plot_dual(time, (flight[:,3] - mars.r)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
AT.plot_dual(time, motion.a_s, flight[:,0], 'Time [s]', 'Acceleration [m/s$^2$]', 'Velocity [m/s]')
AT.plot_single(time, np.degrees(flight[:,2]), 'Time [s]', 'heading angle [deg]')
AT.plot_single(time, np.degrees(flight[:,1]), 'Time [s]', 'flight path angle [deg]')