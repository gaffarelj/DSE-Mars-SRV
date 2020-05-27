import numpy as np
from matplotlib import pyplot as plt
import astro_tools


mean_radius = 3389500                           #[m]
equatorial_radius = 3396200                     #[m]
mu = 42828e9                                    #[m^3/s^2]
scale_height = 8.8e3                            #[m]
rotational_rate = 2*np.pi/(24.6229*3600)        #[rad/s]
J2 = 0.001960454
rho_0 = 0.01417111                              #[kg/m^3]

reentry_altitude = 80000                        #[m]
v_insertion_orbit = 3272.46591                  #[m/s] 
insertion_orbit_a = mean_radius + 500000        #[m]
insertion_orbit_p = mean_radius - 300000        #[m]

cd = 2.75                                       #drag coefficient based on diameter
cl = 1.5                                        #lift coefficicent based on diameter
vehicle_mass = 20000                            #[kg]
S = 65                                          #[m^2]
dt = 0.1 

mars = astro_tools.Planet(mean_radius, scale_height, rho_0, mu, equatorial_radius, J2, rotational_rate)

state = np.zeros(6)
state[0] = mars.reentry_velocity(reentry_altitude, insertion_orbit_a, insertion_orbit_p)                        # velocity
state[1] = mars.reentry_angle(reentry_altitude, insertion_orbit_p, insertion_orbit_a)                           # flight path angle (1.5 arcsecs (2 sigma))
state[2] = np.radians(42.5)                                                                                     # heading angle (1.5 arcsecs (2 sigma))
state[3] = reentry_altitude + mean_radius                                                                       # radius 
state[4] = np.radians(-27.088)                                                                                  # lognitude 25.5
state[5] = np.radians(4.51)                                                                                     # latitude 42.5


motion = astro_tools.Motion(state, np.radians(0), np.radians(0), S, vehicle_mass, cl, cd, mars)
flight, time = motion.forward_euler(dt)
#random = astro_tools.Montecarlo(motion, state, dt)
#random.get_trajectories_linux()

astro_tools.plot_dual(time, (flight[:,3] - mean_radius)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
astro_tools.plot_single(time , np.degrees(flight[:,5]), 'Time [s]', 'latitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,4]), 'Time [s]', 'longitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,2]), 'Time [s]', 'heading angle [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,1]), 'Time [s]', 'flight path angle [deg]')

#astro_tools.scatter(state, random.per)
