import numpy as np
from matplotlib import pyplot as plt
import astro_tools
import csv

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
import aero_calcs as ac

mean_radius = 3389500							#[m]
reentry_altitude = 80000                        #[m]
v_insertion_orbit = 3272.46591                  #[m/s] 
insertion_orbit_a = mean_radius + 500000        #[m]
insertion_orbit_p = mean_radius - 300000        #[m]

vehicle_mass = 38000                            #[kg]
S = 350                                         #[m^2]
MOI =  [391281.1045, 22714482.967, 22714482.967]    # Ixx, Iyy, Izz
dt = 0.1 

mars = astro_tools.Planet()

state = np.zeros(12)
state[0] = mars.reentry_velocity(reentry_altitude, insertion_orbit_a, insertion_orbit_p)                        # velocity
state[1] = mars.reentry_angle(reentry_altitude, insertion_orbit_p, insertion_orbit_a)                           # flight path angle (1.5 arcsecs (2 sigma))
state[2] = np.radians(42.5)                                                                                     # heading angle (1.5 arcsecs (2 sigma))
state[3] = reentry_altitude + mean_radius                                                                       # radius 
state[4] = np.radians(-1.168 - 0.141)                                                                                  # lognitude 25.5 -1.168
state[5] = np.radians(24.322 - 0.077)                                                                                  # latitude 42.5 24.322
state[6] = 0																									# rollrate 
state[7] = 0																									# pitchrate
state[8] = 0																									# yawrate
state[9] = -np.radians(55) 																						# angle of attack defined positive downwards 
state[10] = 0																									# sideslip angle  
state[11] = np.radians(0)																									# bank angle 	


motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars, True, 230e3, 733.45)
flight, time = motion.forward_euler(dt)

astro_tools.plot_dual(time, (flight[:,3] - mean_radius)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
astro_tools.plot_single(time , np.degrees(flight[:,5]), 'Time [s]', 'latitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,4]), 'Time [s]', 'longitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,2]), 'Time [s]', 'heading angle [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,1]), 'Time [s]', 'flight path angle [deg]')
astro_tools.plot_dual(time, motion.mach, motion.q_s, 'Time [s]', 'Mach Number [-]', 'Dynamic Pressure [Pa]')
astro_tools.surfaceplots( np.degrees(flight[:,4]),  np.degrees(flight[:,5]), (flight[:,3] - mean_radius)/1000, 'longitude [deg]', 'latitude [deg]', 'Altitude [km]')

#astro_tools.plot_single(time , -np.degrees(flight[:,9]), 'Time [s]', 'AoA [deg]')
#astro_tools.plot_single(time, motion.roll, 'Time [s]', 'Pitching moment [Nm]')
'''
random = astro_tools.Montecarlo(motion, state, dt)
random.get_trajectories_linux()
astro_tools.scatter(random.per, np.radians(42.5), np.radians(25.5), mars)

with open('accelerations.csv', mode='w') as file:
    accelerations = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    accelerations.writerow(['time', 'downward a', 'forwards a', 'x', 'y', 'z'])
    for i in range(len(motion.a_s)):
        x = (mars.r+reentry_altitude)*np.cos(flight[i,5])*np.cos(flight[i,4])
        y = (mars.r+reentry_altitude)*np.cos(flight[i,5])*np.sin(flight[i,4])
        z = (mars.r+reentry_altitude)*np.sin(flight[i,5])   
        accelerations.writerow([time[i], np.sin(flight[i,1])*motion.a_s[i], np.cos(flight[i,1])*motion.a_s[i], x, y, z])

with open('initial.csv', mode='w') as file:
    initial = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    initial.writerow([ 'downwards a', 'forwards a', 'downwards v', 'forwards v'])
    initial.writerow([np.sin(flight[0,1])*motion.a_s[0], np.cos(flight[0,1])*motion.a_s[0], np.sin(flight[0,1])*flight[0,0], np.cos(flight[0,1])*flight[0,0]])


times = np.arange(740,735,-0.1)
for t in times:
    motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars, True, 230e3, t)
    flight, time = motion.forward_euler(dt)
    if flight[-1,0] < 3:
        print(t)
        break


angles = np.linspace(-5,5,100)
lat = []
long = []

for a in angles:
    state[11] = np.radians(a)
    motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars)          #, True, 230e3, 738.5
    flight, time = motion.forward_euler(dt)
    lat.append(flight[-1,5])
    long.append(flight[-1,4])

plt.scatter(np.degrees(lat), np.degrees(long))
plt.show()

'''