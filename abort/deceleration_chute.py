import sys
sys.path.append(".")
from astrodynamics import mars_standard_atmosphere as MSA
from matplotlib import pyplot as plt

class pc():
	def __init__(self, c_d, A, m):
		self.c_d = c_d
		self.A = A
		self.m = m

class deceleration():
	def __init__(self, V0, h0, chutes, dt=0.1):
		capsule = chutes[0]
		self.m = sum([c.m for c in chutes])
		self.V, self.h, self.t = [V0], [h0], [0]
		self.r = [self.get_rho()]
		self.drag_area = sum([c.c_d * c.A for c in chutes])
		quit = False
		while not quit:# 1: check signs ! velocity starts up (+) and h increases, then V goes down, and so does h
			D = self.comp_drag() - self.m * self.g()
			a_d = D / self.m
			dV = a_d * dt
			dh = (self.V[-1] - dV / 2) * dt
			if self.h[-1] > dh:
				self.V.append(self.V[-1] - dV)
				self.h.append(self.h[-1] - dh)
				self.t.append(self.t[-1] + dt)
				self.r.append(self.rho)
			else:
				quit = True

	def get_rho(self):
		p = MSA.get_pressure(self.h[-1])
		T = MSA.get_temperature(self.h[-1])
		self.rho = MSA.get_density(p, T)
		return self.rho

	def g(self):
		M = 0.64171*10**24
		G = 6.67408*10**(-11)
		R = 3389.5*1000 + self.h[-1]
		return M*G/R**2

	def comp_drag(self):
		return 0.5 * self.get_rho() * self.V[-1] ** 2 * self.drag_area


capsule = pc(0.22, 11.9, 11500)#wrong mass
ballute = pc(0.4, 40, 25)
chutes = [capsule, ballute]
V_0 = 1000
h_0 = 30 * 1e3
# mass, velocity, altitude as function of time
# not only h (x) but also x/y -> longer into atmosphere
decel = deceleration(V_0, h_0, chutes)

plt.rcParams.update({'font.size': 13})
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Time [sec]')
ax1.set_ylabel('Altitude [km]', color=color)
ax1.plot(decel.t, [h/1000 for h in decel.h], color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Velocity [m/s]', color=color)
ax2.plot(decel.t, decel.V, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()
plt.grid()
plt.show()