import numpy as np
import multiprocessing as mp
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import mars_standard_atmosphere as atm
import thermo
import csv
from scipy import interpolate


class Planet:
    def __init__(self, mean_radius=3389500, scale_height=11.1e3, rho_0=0.01417111, gravitational_parameter=42828e9, equatorial_radius=3396200, J2=0.001960454, rotational_rate=7.08824e-5):
        self.r = mean_radius
        self.req = equatorial_radius
        self.mu = gravitational_parameter
        self.J2 = J2
        self.rho_0 = rho_0
        self.omega = rotational_rate
        self.hs = scale_height

    def g(self, altitude, latitude):
        roll_rate = 3 / 2 * np.sin(latitude) ** 2 - 1 / 2
        g = (
            self.mu
            / (self.r + altitude) ** 2
            * (1 - 3 * self.J2 * (self.req / (self.r + altitude)) ** 2 * roll_rate)
        )
        return g

    def density(self, altitude):
        return self.rho_0 * np.exp(-altitude / self.hs)

    def reentry_velocity(self, reentry_altitude, insertion_orbit_a, insertion_orbit_p):
        a_reentry = (insertion_orbit_a + insertion_orbit_p) / 2
        r_reentry = self.r + reentry_altitude
        v_reentry = np.sqrt(self.mu * (2 / r_reentry - 1 / a_reentry))

        return v_reentry

    def reentry_angle(self, reentry_altitude, insertion_orbit_a, insertion_orbit_p):
        a = (insertion_orbit_a + insertion_orbit_p) / 2
        c = a - insertion_orbit_p
        r_reentry = self.r + reentry_altitude
        
        def r(theta):
            return a*(1-(c/a)**2)/(1-(c/a)*np.cos(theta)) - r_reentry

        t = fsolve(r, 0.001)[0]
        x = r_reentry*np.cos(t)
        dydx = - x/(r_reentry*np.sqrt(1 - x**2/r_reentry**2))
        
        return - np.arctan(dydx)

class Motion:
    def __init__(self, inital_conditions, MOI, S, mass, coefficients, Planet, pitch_control = True, thrust = 0, thrust_start = 0, rotation_start = 100000, parachutes=[]):
        self.initial = inital_conditions
        self.Ixx = MOI[0]
        self.Iyy = MOI[1]
        self.Izz = MOI[2]
        self.Planet = Planet
        self.omega = self.Planet.omega
        self.S = S
        self.mass = mass
        self.coefficients = coefficients
        self.chutes = parachutes
        self.i_chute = -1
        self.thrust = thrust
        self.thrust_start = thrust_start
        self.pitch_control = pitch_control
        self.rotation_start = rotation_start

    def dynamicpressure(self, V, r):
        altitude = r - self.Planet.r
        rho = self.Planet.density(altitude)
        return 0.5 * rho * V * V

    def stagnation_heating(self, nose_radius, r,  mach):
        #c = 1.9027e-4
        #q = rho_inf**0.5 * V**3 * c * nose_radius**(-0.5)
        t_wall = 400

        altitude = r - self.Planet.r
        t_static = atm.get_temperature(altitude)
        p_static = atm.get_pressure(altitude)
        if p_static == 0:
            q, t2, rho2 = 0, t_static, 0
        else:
            p2, rho2, t2, m2 = self.normalshock(r, mach)
            gas_wall = thermo.Chemical('carbon dioxide', T=t_wall, P=p2)
            gas_post = thermo.Chemical('carbon dioxide', T=t_wall, P=p2)
            rho_wall = p2/atm.R/t_wall
            mu_wall = gas_wall.mug
            mu_post = gas_post.mug
            t_total = t2*(1+0.5*(atm.gamma - 1)*m2**2)
            dudx = 1/nose_radius * np.sqrt(2*(p2 - p_static)/rho2)
            q = 0.94*(rho2*mu_post)**0.4 * (rho_wall*mu_wall)**0.1 * np.sqrt(dudx) * gas_post.Cpg*(t_total-t_wall)
    
        return q, t2, rho2


    def gravitational_acceleeration(self, r, delta):
        altitude = r - self.Planet.r
        return self.Planet.g(altitude, delta)

    def normalshock(self, r, mach_initial):
        altitude = r - self.Planet.r
        gamma = atm.gamma
        t_static = atm.get_temperature(altitude)
        p_static = atm.get_pressure(altitude)
        rho_static = atm.get_density(p_static, t_static)

        mach = np.sqrt((1+((gamma-1)/2*mach_initial**2))/(gamma*mach_initial**2-(gamma-1)/2))
        pressure = p_static * (1 + 2*gamma/(gamma+1) * (mach_initial**2 - 1))
        density = rho_static * (((gamma+1)*mach_initial**2)/(2+(gamma-1)*mach_initial**2))
        temperature = t_static * (1+(gamma-1)/2*mach_initial**2)*(2*gamma/(gamma-1)*mach_initial**2 - 1)/(mach_initial**2*(2*gamma/(gamma-1) + (gamma-1)/2))

        t_static * pressure/p_static * rho_static/density

        return pressure, density, temperature, mach

    def dVdt(self, g, D, r, gamma, delta, xi):
        dVdt = (
            -D / self.mass
            - g * np.sin(gamma)
            + self.omega ** 2
            * r
            * np.cos(delta)
            * (
                np.sin(gamma) * np.cos(delta)
                - np.cos(gamma) * np.sin(delta) * np.cos(xi)
            )
        )
        return dVdt

    def dgammadt(self, g, D, L, V, r, gamma, delta, xi, mu):
        dgammadt = (
            L * np.cos(mu) / self.mass
            - g * np.cos(gamma)
            + 2 * self.omega * V * np.cos(delta) * np.sin(xi)
            + V * V / r * np.cos(gamma)
            + self.omega ** 2
            * r
            * np.cos(delta)
            * (
                np.cos(gamma) * np.cos(delta)
                - np.sin(gamma) * np.sin(delta) * np.cos(xi)
            )
        ) / V
        return dgammadt

    def dxidt(self, g, L, V, r, gamma, delta, xi, mu):
        dxidt = (
            L * np.sin(mu) / self.mass
            + 2
            * self.omega
            * V
            * (
                np.sin(delta) * np.cos(gamma)
                - np.cos(delta) * np.sin(gamma) * np.cos(xi)
            )
            + V * V / r * np.cos(gamma) ** 2 * np.tan(delta) * np.sin(xi)
            + self.omega ** 2 * r * np.cos(delta) * np.sin(delta) * np.sin(xi)
        ) / (V * np.cos(gamma))
        return dxidt

    def drdt(self, V, gamma):
        return V * np.sin(gamma)

    def dtaudt(self, V, r, gamma, delta, xi):
        return V * np.cos(gamma) * np.sin(xi) / (r * np.cos(delta))

    def ddeltadt(self, V, r, gamma, xi):
        return V / r * np.cos(gamma) * np.cos(xi)

    def dpdt(self, Mx, pitch_rate, yaw_rate):
        return Mx/self.Ixx + (self.Iyy-self.Izz)/self.Ixx * pitch_rate * yaw_rate
    
    def dqdt(self, My, roll_rate, yaw_rate):
        return My/self.Iyy + (self.Izz-self.Ixx)/self.Iyy * roll_rate * yaw_rate
    
    def drratedt(self, Mz, roll_rate, pitch_rate):
        return Mz/self.Izz + (self.Ixx-self.Iyy)/self.Izz * roll_rate * pitch_rate

    def dalphadt(self, roll_rate, pitch_rate, yaw_rate, g, L, V, gamma, mu, alpha, beta):
        return pitch_rate - (roll_rate*np.cos(alpha)+yaw_rate*np.sin(alpha))*np.tan(beta) - (L - self.mass*g*np.cos(gamma)*np.cos(mu))/(self.mass*V*np.cos(beta))
    
    def dbetadt(self, roll_rate, yaw_rate, g, V, gamma, mu, alpha):
        return roll_rate*np.sin(alpha) - yaw_rate*np.cos(alpha) - (self.S+self.mass*g*np.cos(gamma)*np.sin(mu))/(self.mass*V)

    def dmudt(self, roll_rate, q, yaw_rate, g, L, V, gamma, mu, alpha, beta):
        return -(roll_rate*np.cos(alpha) + yaw_rate*np.sin(alpha))/np.cos(beta) - (L - self.mass*g*np.cos(gamma)*np.cos(mu))/(self.mass*V)*np.tan(beta) + (L*np.sin(mu) + self.S*np.cos(mu))/(self.mass*V)*np.tan(gamma)

    def forward_euler(self, timestep):
        flight = [self.initial]
        time = [0]
        self.a_s, self.q_s, self.mach, self.pitch, self.roll, self.heatflux, self.temperature, self.density, self.p2 = [], [], [], [], [], [], [], [], []
        Mx = 0
        My = 0
        Mz = 0
        i = 0
        data = np.genfromtxt('spirit_flightpath.csv', delimiter='', dtype=None)
        ca =abs(data[:,32])
        cn =abs(data[:,34])
        rho = data[:,42]
        aoa = np.radians(data[:,36])
        while flight[-1][3] > self.Planet.r - 3000:
            V          = flight[-1][0]
            gamma      = flight[-1][1]
            xi         = flight[-1][2]
            r          = flight[-1][3]
            tau        = flight[-1][4]
            delta      = flight[-1][5]
            roll_rate  = flight[-1][6]
            pitch_rate = flight[-1][7]
            yaw_rate   = flight[-1][8]
            alpha      = flight[-1][9]
            beta       = flight[-1][10]
            mu         = flight[-1][11]

            #if rho[i] < 0:
                #rho[i] = 0
            q = self.dynamicpressure(V, r)
            altitude = r - self.Planet.r
            mach = V/np.sqrt(atm.gamma*atm.R*atm.get_temperature(altitude))
            postshock_pressure = self.normalshock(r, mach)[0]
            g = self.gravitational_acceleeration(r, delta)

            # Parachute deployment
            chute_drag_area = 0
            if len(self.chutes) > 0:
                # If there's still parachutes after the current one, check if deployment time reached
                if self.i_chute + 1 < len(self.chutes) and time[-1] > self.chutes[self.i_chute+1].deploy_time:
                    # Increment parachute id -> deploy the parachute
                    self.i_chute += 1
                    print(f"Deployed chute {self.i_chute} at {round(time[-1], 2)}")
                    # Remove the mass of the previous parachute from the capsule
                    if self.i_chute >= 1:
                        self.mass -= self.chutes[self.i_chute - 1].m
                # Compute the drag * area of the current parachute
                if self.i_chute >= 0:
                    chute_drag_area = self.chutes[self.i_chute].drag_area

            
            #cl = cn[i]*np.cos(aoa[i]) + ca[i]*np.sin(aoa[i])
            #cd = cn[i]*np.sin(aoa[i]) + ca[i]*np.cos(aoa[i])
            cl,cd = cl_cd(mach, -np.degrees(alpha))

            #if cn[i] == -1:
                #cl = 0
                #cd = 1.8

            D = q * (cd * self.S + chute_drag_area)
            L = q * cl * self.S 
            My = -(L*np.cos(alpha) + D*np.sin(alpha))*0.1
            #i+=1
            if time[-1] > self.thrust_start:
                D = D + self.thrust

            new_state = np.zeros(12)
            a = self.dVdt(g, D, r, gamma, delta, xi)
            new_state[0] = V + timestep * a
            new_state[1] = gamma + timestep * self.dgammadt(g, D, L, V, r, gamma, delta, xi, mu)
            new_state[2] = xi + timestep * self.dxidt(g, L, V, r, gamma, delta, xi, mu)
            new_state[3] = r + timestep * self.drdt(V, gamma)
            new_state[4] = tau + timestep * self.dtaudt(V, r, gamma, delta, xi)
            new_state[5] = delta + timestep * self.ddeltadt(V, r, gamma, xi)
            new_state[6] = roll_rate + timestep * self.dpdt(Mx, pitch_rate, yaw_rate)
            new_state[7] = pitch_rate + timestep * self.dqdt(My, roll_rate, yaw_rate)
            new_state[8] = yaw_rate + timestep * self.drratedt(Mz, roll_rate, pitch_rate)
            new_state[9] = alpha + timestep * self.dalphadt(roll_rate, pitch_rate, yaw_rate, g, L, V, gamma, mu, alpha, beta)
            new_state[10] = beta + timestep * self.dbetadt(roll_rate, yaw_rate, g, V, gamma, mu, alpha)
            new_state[11] = mu + timestep * self.dmudt(roll_rate, q, yaw_rate, g, L, V, gamma, mu, alpha, beta)
            
            if self.pitch_control == True:
                q_dot = (pitch_rate + (roll_rate*np.cos(alpha) + yaw_rate*np.sin(alpha))*np.tan(beta) + (L - self.mass*g*np.cos(gamma)*np.cos(mu))/(self.mass*V*np.cos(beta)) - 
                        new_state[7] - (new_state[6]*np.cos(new_state[9]) + new_state[8]*np.sin(new_state[9]))*np.tan(new_state[10]) - (L - self.mass*g*np.cos(new_state[1])*np.cos(new_state[11]))/(self.mass*V*np.cos(new_state[10]))) / timestep
                p_dot = (-(- yaw_rate*np.sin(alpha) + np.cos(beta)*(- (L - self.mass*g*np.cos(gamma)*np.cos(mu))/(self.mass*V)*np.tan(beta) + (L*np.sin(mu) + self.S*np.cos(mu))/(self.mass*V)*np.tan(gamma)))/np.cos(alpha) + 
                        (- yaw_rate*np.sin(new_state[9]) + np.cos(new_state[10])*(- (L - self.mass*g*np.cos(new_state[1])*np.cos(new_state[11]))/(self.mass*V)*np.tan(new_state[10]) + (L*np.sin(new_state[11]) + self.S*np.cos(new_state[11]))/(self.mass*V)*np.tan(new_state[1])))/np.cos(new_state[9]))
                pitching_moment = q_dot * self.Iyy 
                rolling_moment =  p_dot * self.Izz
                new_state[9] = self.initial[9]
                new_state[11] = self.initial[11]


            self.a_s.append(a)
            self.q_s.append(q)
            self.mach.append(mach)
            self.pitch.append(pitching_moment)
            self.roll.append(rolling_moment)
            q_in, t2, rho2 = self.stagnation_heating(7, r, mach)
            self.heatflux.append(q_in)
            self.temperature.append(t2)
            self.density.append(rho2)
            self.p2.append(postshock_pressure)

            flight.append(new_state)
            time.append(time[-1] + timestep)
            state = new_state

        self.a_s.append(a), self.q_s.append(q), self.mach.append(mach), self.pitch.append(pitching_moment), self.roll.append(rolling_moment), self.heatflux.append(q_in), self.temperature.append(t2), self.density.append(rho2)#, self.p2.append(postshock_pressure)
        return np.array(flight), time

class pc():
    def __init__(self, cd, A, m, deploy_time=0, n=1):
        self.cd = cd
        self.A = A
        self.m = m
        self.n = n
        self.drag_area = A * n * cd
        self.deploy_time = deploy_time

class Montecarlo:
    def __init__(self, Motion, inital_conditions, dt, samples=100):
        self.Motion = Motion
        self.initial = inital_conditions
        self.scale_height = self.Motion.Planet.hs
        self.n = samples
        self.per = None
        self.dt = dt

    def trajectories(self, n):
        np.random.seed(n)
        self.Motion.initial[0] = np.random.normal(self.initial[0], 0.25)                     # 0.25 m/s
        self.Motion.initial[1] = np.random.normal(self.initial[1], np.radians(1 / 60))       # 1 arcsecs
        self.Motion.initial[2] = np.random.normal(self.initial[2], np.radians(1 / 60))       # 1 arcsecs
        self.Motion.initial[3] = np.random.normal(self.initial[3], 8)                        # 8 m
        self.Motion.initial[4] = np.random.normal(self.initial[4], np.radians(7.4313e-4))    
        self.Motion.initial[5] = np.random.normal(self.initial[5], np.radians(5.5487e-3))    
        #self.Motion.initial[9] = np.random.normal(self.initial[9], np.radians(1 / 60))      # 1 arcsecs
        self.Motion.initial[10] = np.random.normal(self.initial[10], np.radians(1 / 60))     # 1 arcsecs
        #self.Motion.initial[11] = np.random.normal(self.initial[11], np.radians(1 / 60))    # 1 arcsecs
        #self.Motion.Planet.hs = np.random.normal(self.scale_height, 50)                     # scale height of atmosphere

        flight, time = self.Motion.forward_euler(self.dt)
        print(n)

        return flight, time

    def impact_point(self, n):
        flight, _ = self.trajectories(n)
        impact = flight[-1]
        return impact

    def get_trajectories_linux(self):
        pool = mp.Pool(mp.cpu_count())
        self.per = pool.map(self.impact_point, range(self.n))


def plot_single(x_data, y_data, x_label, y_label):
    plt.rcParams.update({"font.size": 12})
    fig, ax1 = plt.subplots()

    color = "tab:blue"
    ax1.set_ylabel(y_label)
    ax1.set_xlabel(x_label)
    ax1.plot(x_data, y_data, color=color)
    ax1.tick_params(axis="y")
    #ax1.set_xticks(np.arange(0, 1750, 250))

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.grid()
    plt.show()


def plot_dual(x_data, y_data_1, y_data_2, x_label, y_label_1, y_label_2):
    plt.rcParams.update({"font.size": 12})
    fig, ax1 = plt.subplots()

    color = "tab:red"
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label_1, color=color)
    ax1.plot(x_data, y_data_1, color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    #ax1.set_xticks(np.arange(0, 1750, 250))

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = "tab:blue"
    ax2.set_ylabel(y_label_2, color=color)  # we already handled the x-label with ax1
    ax2.plot(x_data, y_data_2, color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.grid()
    plt.show()

def scatter(stochastic_data, base_lat, base_lon, Planet):
    lat, lon = [], []
    x_distance, y_distance = [], []
    velocity = []

    for i in range(len(stochastic_data)):
        lat.append(stochastic_data[i][5])
        lon.append(stochastic_data[i][4])
        velocity.append(stochastic_data[i][0])

        y = (lat[-1] - base_lat)*Planet.r
        x = (lon[-1] - base_lon)*Planet.r
        x_distance.append(x)
        y_distance.append(y)
    
    with open('impact_points.csv', mode='w') as file:
        impact = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        impact.writerow(['lat', 'long'])
        for i in range(len(lat)):
            impact.writerow([lat[i], lon[i]])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(np.degrees(np.array(lon)), np.degrees(np.array(lat)), s=2)

    confidence_ellipse(np.degrees(np.array(lon)), np.degrees(np.array(lat)), ax, n_std=1,
                   label=r'$1\sigma$', edgecolor='red')
    confidence_ellipse(np.degrees(np.array(lon)), np.degrees(np.array(lat)), ax, n_std=2,
                   label=r'$2\sigma$', edgecolor='fuchsia', linestyle='--')
    confidence_ellipse(np.degrees(np.array(lon)), np.degrees(np.array(lat)), ax, n_std=3,
                   label=r'$3\sigma$', edgecolor='blue', linestyle=':')

    ax.scatter(25.5, 42.5, c='red', s=4)

    ax.axvline(25.5, c='grey', lw=1)
    ax.axhline(42.5, c='grey', lw=1)
    ax.set_title('Impact Points')
    ax.legend()
    plt.ylabel('Lattitude [deg]')
    plt.xlabel('Longitude [deg]')
    plt.savefig('montecarlo_10000')
    plt.show()

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    print(transf)
    return ax.add_patch(ellipse)


def surfaceplots(xdata, ydata, zdata, x_label, y_label, z_label):
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    ax.plot3D(xdata, ydata, zdata, 'red')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    plt.show()


def plot_from_csv(file, x, y):
    data = np.genfromtxt(file, delimiter=",", dtype=None, skip_header = 1)
    lat = data[:,0]
    lon = data[:,1]


    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(np.degrees(np.array(lon)), np.degrees(np.array(lat)), s=2)

    confidence_ellipse(np.degrees(np.array(lon)), np.degrees(np.array(lat)), ax, n_std=1,
                   label=r'$1\sigma$', edgecolor='red')
    confidence_ellipse(np.degrees(np.array(lon)), np.degrees(np.array(lat)), ax, n_std=2,
                   label=r'$2\sigma$', edgecolor='fuchsia', linestyle='--')
    confidence_ellipse(np.degrees(np.array(lon)), np.degrees(np.array(lat)), ax, n_std=3,
                   label=r'$3\sigma$', edgecolor='blue', linestyle=':')

    ax.scatter(25.5, 42.5, c='green', s=4, label = 'mars base')

    ax.scatter(y, x, c='red', s=14, label = 'manoeuvring capability')

    ax.axvline(25.5, c='grey', lw=1)
    ax.axhline(42.5, c='grey', lw=1)
    ax.set_title('Impact Points')
    ax.legend()
    plt.ylabel('Lattitude [deg]')
    plt.xlabel('Longitude [deg]')
    plt.show()


cl_data = np.genfromtxt('cl_standard_config.csv', delimiter=";", dtype=None)
cd_data = np.genfromtxt('cd_standard_config.csv', delimiter=";", dtype=None)
alpha_list = np.round(cl_data[0][1:],1) 
mach_list = np.round(cl_data[:,0],1)
f_cl = interpolate.interp2d(alpha_list, mach_list, cl_data[:, 1:], kind='cubic')
f_cd = interpolate.interp2d(alpha_list, mach_list, cd_data[:, 1:], kind='cubic')

def cl_cd(mach,alpha):
    if mach > 20:
        mach = 20
    if alpha < 40:
        alpha = 40
    if alpha > 60:
        alpha = 60
    if mach < 1.1:
        mach = 1.1
    col = np.where(alpha_list == round(alpha,1))[0][0]+1
    row = np.where(mach_list == round(mach,1))[0][0]

    return f_cl(alpha,mach), f_cd(alpha, mach)#cl_data[row][col], cd_data[row][col]

