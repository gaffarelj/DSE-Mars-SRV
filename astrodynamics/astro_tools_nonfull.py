import numpy as np
import multiprocessing as mp
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

class Planet:
	def __init__(self, mean_radius=3389500, scale_height=11.1e3, rho_0=0.01417111, gravitational_parameter=42828e9, equatorial_radius=3396200, J2=0.001960454, rotational_rate=2 * np.pi / (24.6229 * 3600)):
		self.r = mean_radius
		self.req = equatorial_radius
		self.mu = gravitational_parameter
		self.J2 = J2
		self.rho_0 = rho_0
		self.omega = rotational_rate
		self.hs = scale_height

	def g(self, altitude, latitude):
		p = 3 / 2 * np.sin(latitude) ** 2 - 1 / 2
		g = (self.mu / (self.r + altitude) ** 2 * (1 - 3 * self.J2 * (self.req / (self.r + altitude)) ** 2 * p))
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
			return a * (1 - (c / a) ** 2) / (1 - (c / a) * np.cos(theta)) - r_reentry

		t = fsolve(r, 0.001)[0]
		x = r_reentry * np.cos(t)
		dydx = - x / (r_reentry * np.sqrt(1 - x ** 2 / r_reentry ** 2))
		
		return - np.arctan(dydx)

class Motion:
	def __init__(self, inital_conditions, roll_angle, alpha, S, mass, cl, cd, Planet, parachutes=[], print_deploy=False, prop_reentry=[], end_t=float("Inf")):
		self.initial = inital_conditions
		self.mu = roll_angle
		self.alpha = alpha
		self.Planet = Planet
		self.omega = self.Planet.omega
		self.S = S
		self.mass = mass
		self.cl = cl
		self.cd = cd
		self.chutes = parachutes
		self.i_chute = -1
		self.print_deploy = print_deploy
		self.prop_reentry = prop_reentry
		self.end_time = end_t

	def dynamicpressure(self, V, r):
		altitude = r - self.Planet.r
		rho = self.Planet.density(altitude)
		return 0.5 * rho * V * V

	def gravitational_acceleeration(self, r, delta):
		altitude = r - self.Planet.r
		return self.Planet.g(altitude, delta)

	def dVdt(self, g, D, r, gamma, delta, xi):
		dVdt = (-D / self.mass - g * np.sin(gamma) + self.omega ** 2 * r * np.cos(delta) * (np.sin(gamma) * np.cos(delta) - np.cos(gamma) * np.sin(delta) * np.cos(xi)))
		# dVdt = -D / self.mass - g * np.sin(gamma)
		return dVdt

	def dgammadt(self, g, D, L, V, r, gamma, delta, xi):
		dgammadt = (L * np.cos(self.mu) / self.mass - g * np.cos(gamma) + 2 * self.omega * V * np.cos(delta) * np.sin(xi) + V * V / r * np.cos(gamma) + self.omega ** 2 * r * np.cos(delta) * (np.cos(gamma) * np.cos(delta) - np.sin(gamma) * np.sin(delta) * np.cos(xi))) / V

		# dgammadt = (V / r - g / V) * np.cos(gamma) + (
		#    L * np.cos(self.mu) - self.S * np.sin(self.mu)
		# ) / self.mass / V
		return dgammadt

	def dxidt(self, g, L, V, r, gamma, delta, xi):
		dxidt = (L * np.sin(self.mu) / self.mass + 2 * self.omega * V * (np.sin(delta) * np.cos(gamma) - np.cos(delta) * np.sin(gamma) * np.cos(xi)) + V * V / r * np.cos(gamma) ** 2 * np.tan(delta) * np.sin(xi) + self.omega ** 2 * r * np.cos(delta) * np.sin(delta) * np.sin(xi)) / (V * np.cos(gamma))

		# dxidt = V / r * np.cos(gamma) * np.tan(delta) * np.sin(
		#    xi
		# ) - (L * np.sin(self.mu) - self.S * np.cos(self.mu)) / self.mass / V /
		# np.cos(
		#    gamma
		# )
		return dxidt

	def drdt(self, V, gamma):
		return V * np.sin(gamma)

	def dtaudt(self, V, r, gamma, delta, xi):
		return V * np.cos(gamma) * np.sin(xi) / (r * np.cos(delta))

	def ddeltadt(self, V, r, gamma, xi):
		return V / r * np.cos(gamma) * np.cos(xi)

	def forward_euler(self, timestep):
		flight = [self.initial]
		time = [0]
		self.a_s, self.q_s = [], []
		apogee, entry_burn = False, True
		state = np.zeros(6)
		while flight[-1][3] > self.Planet.r and self.end_time > time[-1]:
			if self.end_time != float("Inf"):
				progress = time[-1]/self.end_time*100
				print(round(progress, 2), "%", sep="", end="\r")
			V = flight[-1][0]
			gamma = flight[-1][1]
			xi = flight[-1][2]
			r = flight[-1][3]
			tau = flight[-1][4]
			delta = flight[-1][5]
			# Entry burn
			if apogee and len(self.prop_reentry) == 2 and q > self.prop_reentry[1] and entry_burn:
				V -= self.prop_reentry[0]
				entry_burn = False

			# Parachute deployment
			chute_drag_area = 0
			if len(self.chutes) > 0:
				# If there's still parachutes after the current one, check if deployment
				# time reached
				if self.i_chute + 1 < len(self.chutes) and time[-1] >= self.chutes[self.i_chute + 1].deploy_time:
					# Increment parachute id -> deploy the parachute
					self.i_chute += 1
					if self.print_deploy:
						print(f"Deployed chute {self.chutes[self.i_chute].name} ({self.i_chute}) at {round(time[-1], 2)}, q={round(self.q_s[-1], 2)} Pa, v={round(V, 2)}")
					# Remove the mass of the previous parachute from the capsule
					if self.i_chute > 0:
					    self.mass -= self.chutes[self.i_chute - 1].m
				# Compute the drag * area of the current parachute
				if self.i_chute >= 0:
					chute_drag_area = self.chutes[self.i_chute].drag_area

			q = self.dynamicpressure(V, r)
			f = 1 if apogee else 1 / 5
			D = q * (self.cd * self.S + chute_drag_area)
			L = q * self.cl * self.S
			g = self.gravitational_acceleeration(r, delta)

			new_state = np.zeros(6)
			a = self.dVdt(g, D, r, gamma, delta, xi)
			new_state[0] = V + timestep * a
			new_state[1] = gamma + timestep * self.dgammadt(g, D, L, V, r, gamma, delta, xi)
			new_state[2] = xi + timestep * self.dxidt(g, L, V, r, gamma, delta, xi)
			new_state[3] = r + timestep * self.drdt(V, gamma)
			new_state[4] = tau + timestep * self.dtaudt(V, r, gamma, delta, xi)
			new_state[5] = delta + timestep * self.ddeltadt(V, r, gamma, xi)
			if not apogee and new_state[3] < state[3]:
				apogee = True

			self.a_s.append(a)
			self.q_s.append(q)
			flight.append(new_state)
			time.append(time[-1] + timestep)
			state = new_state
		if self.i_chute > -1 and self.i_chute == len(self.chutes) - 1:
			self.mass -= self.chutes[self.i_chute].m
		self.a_s.append(a), self.q_s.append(q)
		print()
		return np.array(flight), time

class pc():
	def __init__(self, cd, D, m, deploy_time=0, n=1, name=""):
		self.cd = cd
		self.A = np.pi * D ** 2 / 4
		self.m = m * (n + 1)
		self.n = n
		self.drag_area = self.A * n * cd
		self.deploy_time = deploy_time
		self.name = name

class Montecarlo:
	def __init__(self, Motion, inital_conditions, dt, samples=100):
		self.Motion = Motion
		self.initial = inital_conditions
		self.scale_height = self.Motion.Planet.hs
		self.n = samples
		self.per = None
		self.dt = dt

	def trajectories(self):
		self.Motion.initial[0] = np.random.normal(self.initial[0], 4)                       # 100 m/s
		self.Motion.initial[1] = np.random.normal(self.initial[1], np.radians(1.5 / 60))    # 1.5 arcsecs
		self.Motion.initial[2] = np.random.normal(self.initial[2], np.radians(1.5 / 60))    # 1.5 arcsecs
		self.Motion.initial[3] = np.random.normal(self.initial[3], np.radians(1.5 / 60))    # 1.5 arcsecs
		self.Motion.initial[4] = np.random.normal(self.initial[4], np.radians(1.5 / 60))    # 1.5 arcsecs
		self.Motion.initial[5] = np.random.normal(self.initial[5], np.radians(1.5 / 60))    # 1.5 arcsecs
		
		self.Motion.Planet.hs = np.random.normal(self.scale_height, 100)                    # scale height of atmosphere

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
	plt.rcParams.update({"font.size": 12})
	fig, ax1 = plt.subplots()

	color = "tab:blue"
	ax1.set_xlabel(x_label)
	ax1.set_ylabel(y_label, color=color)
	ax1.plot(x_data, y_data, color=color)
	ax1.tick_params(axis="y", labelcolor=color)
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


def scatter(initial_data, stochastic_data):
	lat = []
	lon = []
	velocity = []

	for i in range(len(stochastic_data)):
		lat.append(stochastic_data[i][5])
		lon.append(stochastic_data[i][4])
		velocity.append(stochastic_data[i][0])

	plt.scatter(np.degrees(np.array(lon)), np.degrees(np.array(lat)))
	plt.ylabel("latitude [deg]")
	plt.xlabel("longitude [deg]")
	plt.grid()
	plt.show()

def angleofattack(V):
	return np.radians(50)