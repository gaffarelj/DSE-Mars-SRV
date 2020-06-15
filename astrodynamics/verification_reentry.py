import numpy as np
from matplotlib import pyplot as plt
import astro_tools
import csv


vehicle_mass = 836                            #[kg]
S = 5.52                                         #[m^2]
MOI = [1, 1, 1]   # Ixx, Iyy, Izz
dt = 0.25 

data = np.genfromtxt('spirit_flightpath.csv', delimiter='', dtype=None)
t = data[:,0]-data[0,0]
alt = data[:,2]
gamma = data[:,1]
lat = data[:,4]
long = data[:,6]
ca = data[:,32]
cn = data[:,34]
vx = data[:,18]
vy = data[:,20]
vz = data[:,22]
v = np.sqrt(vx**2+vy**2+vz**2)
vrel = data[:,24]
aaxial = data[:,8]
aradial = data[:,10]
mars = astro_tools.Planet(scale_height=11.1e3)



state = np.zeros(12)
state[0] = v[0]						                                                     # velocity
state[1] = np.radians(-12)                                                             # flight path angle 
state[2] = np.radians(86.5)                                                              # heading angle 
state[3] = 3522200                                                                       # radius end at 3400300 m
state[4] = np.radians(340.9)                                                             # lognitude 
state[5] = np.radians(-2.9)                                                              # latitude 
state[6] = 0																			 # rollrate 
state[7] = 0																			 # pitchrate
state[8] = 0																			 # yawrate
state[9] = np.radians(0) 																 # angle of attack defined positive downwards 
state[10] = 0																			 # sideslip angle  
state[11] = np.radians(0)																 # bank angle 	


motion = astro_tools.Motion(state, MOI, S, vehicle_mass, 0, mars)     
flight, time = motion.forward_euler(dt)
plt.rcParams.update({"font.size": 12})
plt.plot(t, (alt - mars.r), label = 'flight data')
plt.plot(time, (flight[:,3] - mars.r), label = 'simulation')
plt.legend(loc = 'best')
plt.xlabel('Time [s]')
plt.ylabel('Altitude [m]')
plt.tight_layout()
plt.show()

plt.plot(t, v, label = 'flight data')
plt.plot(time, flight[:,0], label = 'simulation')
plt.legend(loc = 'best')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.show()

plt.plot(t, lat, label = 'flight data')
plt.plot(time, np.degrees(flight[:,5]), label = 'simulation')
plt.legend(loc = 'best')
plt.xlabel('Time [s]')
plt.ylabel('Latitude [deg]')
plt.show()

plt.plot(t, long, label = 'flight data')
plt.plot(time, np.degrees(flight[:,4]), label = 'simulation')
plt.legend(loc = 'best')
plt.xlabel('Time [s]')
plt.ylabel('Longitude [deg]')
plt.show()

height1 = alt - mars.r
height2 = flight[:,3] - mars.r

e1 = np.dot((height1 - height2),(height1 - height2))/np.dot(height1,height1)
e2 = np.dot((v - flight[:,0]),(v - flight[:,0]))/np.dot(v,v)
e3 = np.dot((lat - np.degrees(flight[:,5])),(lat - np.degrees(flight[:,5])))/np.dot(lat,lat)
e4 = np.dot((long - np.degrees(flight[:,4])),(long - np.degrees(flight[:,4])))/np.dot(long,long)

print(e1*100)
print(e2*100)
print(e3*100)
print(e4*100)

astro_tools.plot_single(time , motion.q_s, 'Time [s]', 'Stagnation heat flux [W/m^2]')

'''

astro_tools.plot_dual(time, (flight[:,3] - mars.r)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
astro_tools.plot_single(time , np.degrees(flight[:,5]), 'Time [s]', 'latitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,4]), 'Time [s]', 'longitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,2]), 'Time [s]', 'heading angle [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,1]), 'Time [s]', 'flight path angle [deg]')
astro_tools.surfaceplots( np.degrees(flight[:,4]),  np.degrees(flight[:,5]), (flight[:,3] - mars.r)/1000, 'longitude [deg]', 'latitude [deg]', 'Altitude [km]' )
'''


