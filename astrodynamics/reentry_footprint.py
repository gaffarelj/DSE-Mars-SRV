import numpy as np
from matplotlib import pyplot as plt
import astro_tools
import mars_standard_atmosphere as atm

mean_radius = 3389500                           #[m]
equatorial_radius = 3396200                     #[m]
mu = 42828e9                                  #[m^3/s^2]
reentry_altitude = 80000                        #[m]
v_reentry = 3523.47726                          #[m/s]
scale_height = 8.8e3                        #[m]
rotational_rate = 2*np.pi/(24.6229*3600)          #[rad/s]
J2 = 0.001960454
rho_0 = 0.01417111                              #[kg/m^3]

cd = 2.75                                       #drag coefficient based on diameter
cl = 1.5                                        #lift coefficicent based on diameter
vehicle_mass = 20000                            #[kg]
S = 65     
dt = 0.1 

mars = astro_tools.Planet(mean_radius, scale_height, rho_0, mu, equatorial_radius, J2, rotational_rate)

state = np.zeros(6)
state[0] = v_reentry                                    # velocity
state[1] = np.radians(-1.5)                             # flight path angle 
state[2] = np.radians(10)                               # heading angle 
state[3] = reentry_altitude + mean_radius               # radius 
state[4] = np.radians(10)                               # lognitude 
state[5] = np.radians(40)                               # latitude

flight = [state]
time = [0]

while flight[-1][3] > mean_radius:
    time.append(time[-1] + dt)
    translation = astro_tools.Motion(state, 0, S, vehicle_mass, cl, cd, mars)
    new_state = translation.forward_euler(dt)
    flight.append(new_state)
    state = new_state

flight = np.array(flight)
astro_tools.plot_dual(time, (flight[:,3] - mean_radius)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
astro_tools.plot_single(time , np.degrees(flight[:,5]), 'Time [s]', 'latitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,4]), 'Time [s]', 'longitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,2]), 'Time [s]', 'heading angle [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,1]), 'Time [s]', 'flight path angle [deg]')
