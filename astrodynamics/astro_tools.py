import numpy as np
import multiprocessing as mp
from matplotlib import pyplot as plt


class Planet:
    def __init__(
        self,
        mean_radius,
        scale_height,
        rho_0,
        gravitational_parameter,
        equatorial_radius,
        J2,
        rotational_rate,
    ):
        self.r = mean_radius
        self.req = equatorial_radius
        self.mu = gravitational_parameter
        self.J2 = J2
        self.rho_0 = rho_0
        self.omega = rotational_rate
        self.hs = scale_height

    def g(self, altitude, latitude):
        p = 3 / 2 * np.sin(latitude) ** 2 - 1 / 2
        g = (
            self.mu
            / (self.r + altitude) ** 2
            * (1 - 3 * self.J2 * (self.req / (self.r + altitude)) ** 2 * p)
        )
        return g

    def density(self, altitude):
        return self.rho_0 * np.exp(-altitude / self.hs)

    def reentry_velocity(self, reentry_altitude, insertion_orbit_a, insertion_orbit_p):
        a_reentry = (insertion_orbit_a + insertion_orbit_p) / 2
        r_reentry = (self.r + reentry_altitude) 
        v_reentry = np.sqrt(self.mu*(2/r_reentry - 1/a_reentry))

        return v_reentry


class Motion:
    def __init__(self, inital_conditions, roll_angle, S, mass, cl, cd, Planet): 
        self.initial = inital_conditions
        self.mu = roll_angle
        self.Planet = Planet
        self.S = S
        self.mass = mass
        self.cl = cl
        self.cd = cd

    def dynamicpressure(self, V, r):
        altitude = r - self.Planet.r
        rho = self.Planet.density(altitude)
        return 0.5 * rho * V * V

    def gravitational_acceleeration(self, r, delta):
        altitude = r - self.Planet.r
        return self.Planet.g(altitude, delta)

    def dVdt(self, g, D, gamma):
        dVdt = -D / self.mass - g * np.sin(gamma)
        return dVdt

    def dgammadt(self, g, D, L, V, r, gamma):
        dgammadt = (V / r - g / V) * np.cos(gamma) + (
            L * np.cos(self.mu) - self.S * np.sin(self.mu)
        ) / self.mass / V
        return dgammadt

    def dxidt(self, g, L, V, r, gamma, delta, xi):
        dxidt = V / r * np.cos(gamma) * np.tan(delta) * np.sin(
            xi
        ) - (L * np.sin(self.mu) - self.S * np.cos(self.mu)) / self.mass / V / np.cos(
            gamma
        )
        return dxidt

    def drdt(self, V, gamma):
        return V * np.sin(gamma)

    def dtaudt(self, V, r, gamma, delta, xi):
        return (
            V
            * np.cos(gamma)
            * np.sin(xi)
            / (r * np.cos(delta))
        )

    def ddeltadt(self, V, r, gamma, xi):
        return V / r * np.cos(gamma) * np.cos(xi)

    def forward_euler(self, timestep):
        flight = [self.initial]
        time = [0]
        while flight[-1][3] > self.Planet.r:

            V = flight[-1][0]
            gamma = flight[-1][1]
            xi = flight[-1][2]
            r = flight[-1][3]
            tau = flight[-1][4]
            delta = flight[-1][5]

            D = self.dynamicpressure(V, r) * self.cd * self.S
            L = self.dynamicpressure(V, r) * self.cl * self.S
            g = self.gravitational_acceleeration(r, delta)

            new_state = np.zeros(6)
            new_state[0] = V + timestep * self.dVdt(g,D,gamma)
            new_state[1] = gamma + timestep * self.dgammadt(g, D, L, V, r, gamma)
            new_state[2] = xi + timestep * self.dxidt(g, L, V, r, gamma, delta, xi)
            new_state[3] = r + timestep * self.drdt(V, gamma)
            new_state[4] = tau + timestep * self.dtaudt(V,r,gamma,delta,xi)
            new_state[5] = delta + timestep * self.ddeltadt(V,r,gamma,xi)

            flight.append(new_state)
            time.append(time[-1] + timestep)
            state = new_state
        return np.array(flight), time
        

class Montecarlo():
    def __init__(self, Motion, inital_conditions, dt, samples = 100):
        self.Motion = Motion
        self.initial = inital_conditions
        self.scale_height = self.Motion.Planet.hs
        self.n = samples
        self.per = None
        self.dt = dt

    def trajectories(self):
        self.Motion.initial[0] = np.random.normal(self.initial[0], 4)                                # 100 m/s
        self.Motion.initial[1] = np.random.normal(self.initial[1], np.radians(1.5/60))               # 1.5 arcsecs
        self.Motion.initial[2] = np.random.normal(self.initial[2], np.radians(1.5/60))               # 1.5 arcsecs
        self.Motion.initial[3] = np.random.normal(self.initial[3], np.radians(1.5/60))               # 1.5 arcsecs
        self.Motion.initial[4] = np.random.normal(self.initial[4], np.radians(1.5/60))               # 1.5 arcsecs
        self.Motion.initial[5] = np.random.normal(self.initial[5], np.radians(1.5/60))               # 1.5 arcsecs

        self.Motion.Planet.hs = np.random.normal(self.scale_height, 100)                             # scale height of atmosphere 

        flight, time = self.Motion.forward_euler(self.dt)

        return flight, time 

    def impact_point(self, n):
        flight, time = self.trajectories()
        impact = flight[-1]
        return impact 

    def get_trajectories_linux(self):
        pool = mp.Pool(mp.cpu_count())
        self.per = pool.map(self.impact_point, range(self.n))



def plot_single(x_data, y_data, x_label, y_label):
    plt.rcParams.update({'font.size': 12})
    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label, color=color)
    ax1.plot(x_data, y_data, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_xticks(np.arange(0, 1750, 250))
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.grid()
    plt.show()
    

def plot_dual(x_data, y_data_1, y_data_2, x_label, y_label_1, y_label_2):
    plt.rcParams.update({'font.size': 12})
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label_1, color=color)
    ax1.plot(x_data, y_data_1, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_xticks(np.arange(0, 1750, 250))

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(y_label_2, color=color)  # we already handled the x-label with ax1
    ax2.plot(x_data, y_data_2, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.grid()
    plt.show()

def scatter(initial_data, stochastic_data):
    lat = []
    lon = []
    velocity = []
    
    for i in range(len(stochastic_data)):
        lat.append(stochastic_data[i][5])
        lon.append(stochastic_data[i][4])
        velocity.append(stochastic_data[i][0])
        
    plt.scatter(np.degrees(np.array(lon)), np.degrees(np.array(lat)))
    plt.ylabel('latitude [deg]')
    plt.xlabel('longitude [deg]')
    plt.grid()
    plt.show()
