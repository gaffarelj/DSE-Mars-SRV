import numpy as np
from matplotlib import pyplot as plt
import mars_standard_atmosphere as atm 

# Constants  
mean_radius = 3389500                           #[m]
mu = 42828                                      #[km^3/s^2]
v_insertion_orbit = 3272.46591                  #[m/s] 
reentry_altitude = 80000                        #[m]
insertion_orbit_a = mean_radius + 500000        #[m]
insertion_orbit_p = mean_radius - 300000        #[m]
gamma = atm.gamma
R = atm.R
g0 = 3.71                                       #[m/s^2]
ct = 2.75                                       #tangential force coefficient based on diameter
cn = 1.5                                        #normal force coefficicent based on diameter
vehicle_mass = 20000                            #[kg]
S = 65                                          #[m^2] 


def reentry_orbit(reentry_altitude, insertion_orbit_a, insertion_orbit_p, mu=mu):
    a_reentry = (insertion_orbit_a + insertion_orbit_p)/2000
    r_reentry = (mean_radius+reentry_altitude)/1000
    v_reentry = np.sqrt(mu*(2/r_reentry - 1/a_reentry))*1000

    return v_reentry

# Dynamic Pressure
def dynamic_pressure(V, altitude, gamma=gamma):
    pressure = atm.get_pressure(altitude)
    temperature = atm.get_temperature(altitude)
    density = atm.get_density(pressure, temperature)
    a= atm.get_speed_of_sound(temperature)
    mach = V/a

    if altitude > 100000:                       #Altitude too high to cause significant dynamic pressure 
        q = 0
    else:
        q = mach**2 * 0.5*gamma*pressure

    return q

def normal_shock(V, pressure_1, density_1, temperature_1):
    a_1 = atm.get_speed_of_sound(temperature_1)
    mach_1 = V/a_1   

    mach_2 = np.sqrt((1+((gamma-1)/2*mach_1**2))/(gamma*mach_1**2-(gamma-1)/2))
    pressure_2 = pressure_1 * (1 + 2*gamma/(gamma+1) * (mach_1**2 - 1))
    density_2 = density_1 * (((gamma+1)*mach_1**2)/(2+(gamma-1)*mach_1**2))
    temperature_2 = temperature_1 * pressure_2/pressure_1 * density_1/density_2

    return pressure_2, density_2, temperature_2, mach_2

# Ballistic equations
def ballistic_coefficient(g, S=S, ct=ct, m=vehicle_mass):
    return m*g/(ct*S)

def gravitational_acceleration(h, mean_radius=mean_radius, g0=g0):
    return g0*(mean_radius/(mean_radius + h))**2

def dVdt(g, q, beta, flight_path_angle):
    return g*(-q/beta + np.sin(flight_path_angle))

def dalphadt(g, q, beta, flight_path_angle, V, h, cn=cn, ct=ct):
    return (-q*g/beta * cn/ct + np.cos(flight_path_angle) * (g-V**2/(mean_radius+h)))/V

def dhdt(V, flight_path_angle):
    return -V*np.sin(flight_path_angle)

def drdt(flight_path_angle, V, h):
    return mean_radius*V*np.cos(flight_path_angle)/(mean_radius+h)

# Integration scheme 
def forward_euler(y_n, dy):
    return y_n + dy


# Initial Conditions 
v_reentry = reentry_orbit(reentry_altitude, insertion_orbit_a, insertion_orbit_p)
height = [reentry_altitude]
velocity = [v_reentry]
flight_path = [1*np.pi/180]
time = [0]
distance = [0]
deceleration = [0]
dyn_pressure = [0]
dt = 0.01

# Reentry trajectory
while height[-1] > 1000:
    time.append(time[-1] + dt)
    
    g = gravitational_acceleration(height[-1])
    beta = ballistic_coefficient(g)
    q = dynamic_pressure(velocity[-1], height[-1])
    dyn_pressure.append(q)

    dV = dVdt(g, q, beta, flight_path[-1]) * dt
    deceleration.append(dV/dt)
    V = forward_euler(velocity[-1], dV)
    velocity.append(V)

    dh = dhdt(velocity[-1], flight_path[-1]) * dt
    h = forward_euler(height[-1], dh)
    height.append(h)

    dalpha = dalphadt(g, q, beta, flight_path[-1], velocity[-1], height[-1]) * dt
    alpha = forward_euler(flight_path[-1], dalpha)
    flight_path.append(alpha)

    dr = drdt(flight_path[-1], velocity[-1], height[-1]) * dt
    r = forward_euler(distance[-1], dr)
    distance.append(r) 


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Altitude [km]', color=color)
ax1.plot(time, np.array(height)/1000, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Velocity [m/s]', color=color)  # we already handled the x-label with ax1
ax2.plot(time, velocity, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.show()
