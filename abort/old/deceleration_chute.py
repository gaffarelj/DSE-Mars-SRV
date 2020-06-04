import sys
sys.path.append(".")
from astrodynamics import mars_standard_atmosphere as MSA
from matplotlib import pyplot as plt
import numpy as np

class pc():
	def __init__(self, c_d, A, m, asc_cd=0):
		self.c_d = c_d
		self.A = A
		self.m = m
		self.asc_cd = asc_cd

	def get_cd(self, V=0):
		if self.asc_cd == 0:
			return self.c_d
		elif V > 0:
			return self.asc_cd
		return self.c_d

class deceleration():
	def __init__(self, V0, h0, chutes, dt=0.1):
		self.chutes = chutes
		self.m = sum([c.m for c in self.chutes])
		self.V, self.h, self.t = [V0], [h0], [0]
		self.r = [self.get_rho()]
		self.D, self.W, self.F, self.a = [], [], [], []
		self.i_chute = 1
		self.comp_drag_area()
		while True:
			D = - np.sign(self.V[-1]) * self.comp_drag()
			W = self.m * self.g()
			F = D - W
			a = F / self.m
			self.D.append(D), self.W.append(W), self.F.append(F), self.a.append(a)
			dV = a * dt
			dh = (self.last_V() + dV / 2) * dt
			if self.last_h() > 0:
				self.V.append(self.last_V() + dV)
				self.h.append(self.last_h() + dh)
				self.t.append(self.t[-1] + dt)
				self.r.append(self.rho)
				if dh < 0 and self.i_chute == 1:
					self.i_chute = 2
					self.comp_drag_area()
			else:
				break

	def comp_drag_area(self):
		self.drag_area = sum([c.get_cd(self.last_V()) * c.A for c in self.chutes[0:self.i_chute]])

	def get_rho(self):
		p = MSA.get_pressure(self.last_h())
		T = MSA.get_temperature(self.last_h())
		self.rho = MSA.get_density(p, T)
		return self.rho

	def g(self):
		M = 0.64171 * 10 ** 24
		G = 6.67408 * 10 ** (-11)
		R = 3389500 + self.last_h()
		return M * G / R ** 2

	def comp_drag(self):
		return 0.5 * self.get_rho() * self.last_V() ** 2 * self.drag_area
	
	def last_V(self):
		return self.V[-1]
	def last_h(self):
		return self.h[-1]


capsule = pc(1.2, 11.9, 11500, asc_cd=0.4)
ballute = pc(0.4, 63.6, 25)
parachute = pc(0.8, 3 * 600, 50)
chutes = [capsule, parachute]
V_0 = 400
h_0 = 10 * 1e3
# 3: mass, velocity, altitude as function of time (for starting point)
# 2: not only h (x) but also x/y -> longer into atmosphere ?
decel = deceleration(V_0, h_0, chutes)

t = 100
t_i = np.where(np.array(decel.t) >= t)[0][0]
plt.rcParams.update({'font.size': 13})

time = decel.t[t_i:]
fig, axs = plt.subplots(2, 2)

axs[0, 0].plot(time, decel.V[t_i:])
axs[0, 0].set_ylabel('Velocity [m/s]')

axs[0, 1].plot(time, [h/1000 for h in decel.h[t_i:]])
axs[0, 1].set_ylabel('Altitude [km]')

axs[1, 0].plot(time, decel.a[t_i:])
axs[1, 0].set_ylabel('Acceleration [m/s$^2$]')

axs[1, 1].plot(time, decel.r[t_i:])
axs[1, 1].set_ylabel('ρ [kg/m$^3$]')

for t in axs:
	for ax in t:
		ax.grid()
		ax.set_xlabel('Time [sec]')

plt.show()