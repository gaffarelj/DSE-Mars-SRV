import numpy as np
from matplotlib import pyplot as plt
import sys
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
import aero_calcs as ac
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'astrodynamics'))
import astro_tools

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
state[4] = np.radians(-1.168)                                                                                  # lognitude 25.5 -1.168 - 0.141
state[5] = np.radians(24.322)                                                                                  # latitude 42.5 24.322 - 0.077
state[6] = 0																									# rollrate 
state[7] = 0																									# pitchrate
state[8] = 0																									# yawrate
state[9] = -np.radians(60) 																						# angle of attack defined positive downwards 
state[10] = 0																									# sideslip angle  
state[11] = np.radians(0)																									# bank angle 	


motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars) #, True, 230e3, 733.45)
flight, time = motion.forward_euler(dt)


astro_tools.plot_single(time , motion.heatflux, 'Time [s]', 'Heat Flux [W/m^2]')