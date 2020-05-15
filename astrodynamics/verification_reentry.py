import numpy as np
from matplotlib import pyplot as plt
import isa

# Constants  
mean_radius = 6731e3                #[m]
reentry_altitude = 400000*0.3048    #[m]
v_reentry = 36309*0.3048            #[m/s]
gamma = 1.4
g0 = 9.81                           #[m/s^2]
R = 287
ct = 1.28                            # tangential force coefficient 
cn = 0.22                           # normal force coefficicent 
vehicle_mass = 13000*0.453592                 #[kg]
S = 3.9**2*np.pi/4             #[m^2] 


def dynamic_pressure_earth(V, altitude, gamma=gamma, R=R):
    temperature_1, pressure_1, density_1 = isa.ISA(altitude)
    a_1 = np.sqrt(gamma*R*temperature_1)
    mach_1 = V/a_1

    q = mach_1**2 * 0.5*gamma*pressure_1

    return q

# Ballistic equations
def ballistic_coefficient(g, S=S, ct=ct, m=vehicle_mass):
    return m*g/ct/S

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


# Reentry trajectory
height = [reentry_altitude]
velocity = [v_reentry]
flight_path = [6.62*np.pi/180]
time = [0]
distance = [0]
deceleration = [0]
dyn_pressure = [0]
dt = 0.01

while height[-1] > 11000:
    time.append(time[-1] + dt)
    
    g = gravitational_acceleration(height[-1])
    beta = ballistic_coefficient(g)
    q = dynamic_pressure_earth(velocity[-1], height[-1])
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

plt.rcParams.update({'font.size': 14})
plt.plot(np.array(time)/60, np.array(height)/1000/0.3048)
plt.ylabel("Altitude [1000 ft]")
plt.xlabel("Time [min]")
plt.grid()
plt.show()

plt.plot(np.array(time)/60, np.array(velocity)/1000/0.3048)
plt.ylabel("Velocity [1000 ft/s]")
plt.xlabel("Time [min]")
plt.grid()
plt.show()


plt.plot(np.array(time)/60, np.array(flight_path)/1000/0.3048)
plt.ylabel("Velocity [1000 ft/s]")
plt.xlabel("Time [min]")
plt.grid()
plt.show()

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Altitude [1000 ft]', color=color)
ax1.plot(np.array(time)/60, np.array(height)/1000/0.3048, color=color)
ax1.set_xticks(np.arange(0, 9, 1))
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Velocity [1000 ft/s]', color=color)  # we already handled the x-label with ax1
ax2.plot(np.array(time)/60, np.array(velocity)/1000/0.3048, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.show()
