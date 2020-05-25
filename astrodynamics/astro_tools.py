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

    def entry_orbit(self):
        return 3


class Motion:
    def __init__(self, inital_conditions, roll_angle, S, mass, cl, cd, Planet): 
        self.initial = inital_conditions
        self.V = inital_conditions[0]
        self.gamma = inital_conditions[1]
        self.xi = inital_conditions[2]
        self.r = inital_conditions[3]
        self.tau = inital_conditions[4]
        self.delta = inital_conditions[5]
        self.mu = roll_angle
        self.Planet = Planet
        self.S = S
        self.mass = mass
        self.cl = cl
        self.cd = cd

    def dynamicpressure(self):
        altitude = self.r - self.Planet.r
        rho = self.Planet.density(altitude)
        return 0.5 * rho * self.V * self.V

    def gravitational_acceleeration(self):
        altitude = self.r - self.Planet.r
        return self.Planet.g(altitude, self.delta)

    def dVdt(self, g, D):
        dVdt = -D / self.mass - g * np.sin(self.gamma)
        return dVdt

    def dgammadt(self, g, D, L):
        dgammadt = (self.V / self.r - g / self.V) * np.cos(self.gamma) + (
            L * np.cos(self.mu) - self.S * np.sin(self.mu)
        ) / self.mass / self.V
        return dgammadt

    def dxidt(self, g, L):
        dxidt = self.V / self.r * np.cos(self.gamma) * np.tan(self.delta) * np.sin(
            self.xi
        ) - (L * np.sin(self.mu) - self.S * np.cos(self.mu)) / self.mass / self.V / np.cos(
            self.gamma
        )
        return dxidt

    def drdt(self):
        return self.V * np.sin(self.gamma)

    def dtaudt(self):
        return (
            self.V
            * np.cos(self.gamma)
            * np.sin(self.xi)
            / (self.r * np.cos(self.delta))
        )

    def ddeltadt(self):
        return self.V / self.r * np.cos(self.gamma) * np.cos(self.xi)

    def forward_euler(self, timestep):
        altitude = self.r - self.Planet.r
        D = self.dynamicpressure() * self.cd * self.S
        L = self.dynamicpressure() * self.cl * self.S
        g = self.gravitational_acceleeration()

        new_state = np.zeros(6)
        new_state[0] = self.V + timestep * self.dVdt(g,D)
        new_state[1] = self.gamma + timestep * self.dgammadt(g, D, L)
        new_state[2] = self.xi + timestep * self.dxidt(g, L)
        new_state[3] = self.r + timestep * self.drdt()
        new_state[4] = self.tau + timestep * self.dtaudt()
        new_state[5] = self.delta + timestep * self.ddeltadt()
        return new_state


class Montecarlo():
    def __init__(self, Iterator, initial_state, parameter, mean, standard_deviation, cutoff_radius, samples = 100):
        self.Iterator = Iterator
        self.initial_state = initial_state
        self.parameter = parameter
        self.mean = mean
        self.sigma = standard_deviation
        self.stop = cutoff_radius
        self.n = samples
        self.per = None

    #def trajectories(self):

    def get_trajectories_linux(self):
        pool = mp.Pool(mp.cpu_count())
        self.per = pool.map(self.trajectories, range(self.n))
        self.per = np.sum(self.per,axis=0)/self.n


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

