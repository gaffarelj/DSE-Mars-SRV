import numpy as np
from matplotlib import pyplot as plt
import astro_tools
import csv
import matplotlib
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
import aero_calcs as ac
import mars_standard_atmosphere as atm

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
state[9] = -np.radians(55) 																						# angle of attack defined positive downwards 
state[10] = 0																									# sideslip angle  
state[11] = np.radians(0)																									# bank angle 	


motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars) #, True, 230e3, 733.45)
flight, time = motion.forward_euler(dt)

astro_tools.plot_dual(time, (flight[:,3] - mean_radius)/1000, flight[:,0], 'Time [s]', 'Altitude [km]', 'Velocity [m/s]')
astro_tools.plot_single(time , np.degrees(flight[:,5]), 'Time [s]', 'latitude [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,4]), 'Time [s]', 'longitude [deg]')

astro_tools.plot_single(time , np.degrees(flight[:,2]), 'Time [s]', 'heading angle [deg]')
astro_tools.plot_single(time , np.degrees(flight[:,1]), 'Time [s]', 'flight path angle [deg]')
astro_tools.plot_dual(time, motion.mach, motion.q_s, 'Time [s]', 'Mach Number [-]', 'Dynamic Pressure [Pa]')
astro_tools.surfaceplots( np.degrees(flight[:,4]),  np.degrees(flight[:,5]), (flight[:,3] - mean_radius)/1000, 'longitude [deg]', 'latitude [deg]', 'Altitude [km]')
astro_tools.plot_single(time , motion.heatflux, 'Time [s]', 'Heat Flux [W/m^2]')
astro_tools.plot_single(time , motion.temperature, 'Time [s]', 'Post Shock Temperature [K]')
astro_tools.plot_single(time , motion.density, 'Time [s]', 'Post Shock Density [kg/m^3]')

#astro_tools.plot_single(time , -np.degrees(flight[:,9]), 'Time [s]', 'AoA [deg]')
astro_tools.plot_single(time, motion.pitch, 'Time [s]', 'Rolling moment [Nm]')

#theta = motion.pitch/np.degrees((19.28-1.7 - 6.89 + 1.7/2)*np.sqrt(pow(motion.mach,2)-1)/7.755/2)
#astro_tools.plot_single(time, theta, 'Time [s]', 'Pitching moment [Nm]')

idx = 7334
m = motion.mach[idx]
q = motion.q_s[idx]
p2 = motion.p2[idx]
alt = flight[idx,3] - mean_radius
rho = mars.density(alt)
p = atm.get_pressure(alt)
q = motion.q_s[idx]
print(motion.pitch[idx],flight[idx,0], rho)
dp = abs(motion.pitch[idx] / (19.28 - 1.7 - 6.89 + 1.7/2))
cp = dp/(q*5.185)


print(dp, cp)

'''
random = astro_tools.Montecarlo(motion, state, dt, 10000)
random.get_trajectories_linux()
astro_tools.scatter(random.per, np.radians(42.5), np.radians(25.5), mars)

max_pitch = np.radians(5)
max_roll = np.radians(3.4)
x, y = [], []

state[9] = -np.radians(55) + max_pitch
motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars)
flight, time = motion.forward_euler(dt)
x.append(np.degrees(flight[-1,5]))# + 0.9266)
y.append(np.degrees(flight[-1,4]))# + 2.2395)


state[9] = -np.radians(55) - max_pitch
motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars)
flight, time = motion.forward_euler(dt)
x.append(np.degrees(flight[-1,5]))
y.append(np.degrees(flight[-1,4]))

state[9] = -np.radians(55)

state[11] = max_roll
motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars)
flight, time = motion.forward_euler(dt)
x.append(np.degrees(flight[-1,5]))
y.append(np.degrees(flight[-1,4]))

state[11] = -max_roll
motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars)
flight, time = motion.forward_euler(dt)
x.append(np.degrees(flight[-1,5]))
y.append(np.degrees(flight[-1,4]))

astro_tools.plot_from_csv('impact_points.csv', x, y)
'''

'''
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
'''
r = mars.r/1000
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0.0:pi/2:50j, 0.0:pi/2:50j]
xp = r*sin(phi)*cos(theta)
yp = r*sin(phi)*sin(theta)
zp = r*cos(phi)


# convert to 2d matrices
Z = np.outer(zp.T, zp)        # 50x50
X, Y = np.meshgrid(xp, yp)    # 50x50

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
X, Y, Z = xp, yp, zp
C = mars.g(0, phi).reshape(Z.shape)
scamap = plt.cm.ScalarMappable(cmap='inferno')
fcolors = scamap.to_rgba(C)
ax.plot_surface(X, Y, Z, facecolors=fcolors, cmap='inferno', label = 'g [m/s^2]')
fig.colorbar(scamap)
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')
plt.show()

'''

