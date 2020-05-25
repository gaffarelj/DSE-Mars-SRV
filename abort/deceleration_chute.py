import sys
sys.path.append(".")
from astrodynamics import mars_standard_atmosphere as MSA
import numpy as np

class pc():
	def __init__(self, c_d, A, m):
		self.c_d = c_d
		self.A = A
		self.m = m

class deceleration():
	def __init__(self, V0, h0, chutes, g=3.711):
		capsule = chutes[0]
		self.m = sum([c.m for c in chutes])
		self.V, self.h, self.t = [V0], [h0], [0]
		self.g = g
		self.get_rho()
		self.drag_area = sum([c.c_d*c.A for c in chutes])

	def get_rho(self):	# in meters
		p = MSA.get_pressure(self.h[-1])
		T = MSA.get_temperature(self.h[-1])
		self.rho = MSA.get_density(p, T)

	def terminal_velocity(self):
		V_ter = np.sqrt(2*self.m*self.g / (self.rho * self.drag_area))
		self.V.append(V_ter)
		self.dV = self.V[-2] - self.V[-1]

	def dec_time(self):
		t = (1/self.V[-1] - 1/self.V[-2])*2*self.m/(self.rho*self.drag_area)
		self.t.append(t)

	def new_h(self):
		print(self.t[-1])
		h = self.h[-1] - 0.5*self.dV*self.t[-1]
		self.h.append(h)

capsule = pc(0.22, 11.9, 11500)
ballute = pc(0.4, 15, 25)
chutes = [capsule, ballute]
V_0 = 4188.75
h_0 = 20 * 1e3

D = deceleration(V_0, h_0, chutes)
V_1 = D.terminal_velocity()
D.dec_time()
D.new_h()
print([t/60 for t in D.t], D.dV, D.V, D.h)