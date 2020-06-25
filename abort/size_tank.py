import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT

class tank():
	def __init__(self, r0, r1, V, P, stress_max=276e6, rho=2700):
		# Min/max radius, in m
		# Required volume
		# Tank pressure, in Pa
		self.r0 = r0
		self.r1 = r1
		self.rho = rho
		self.V_req = V
		self.P, self.s_max = P, stress_max
		self.r, self.L, self.t, self.m = [], [], [], []
		self.run()

	def run(self, L_r_min=3, L_r_max=20):
		r_s = np.arange(self.r0, self.r1, 0.0001)
		L_s = [self.L_r(r) for r in r_s]
		L_neg = np.where(np.array(L_s) <= 0)[0][0]
		r_s, L_s = r_s[:L_neg], L_s[:L_neg]
		L_r_max_i = -np.where(np.flip(np.array(L_s)/np.array(r_s)) >= L_r_max*1.01)[0][0]
		L_r_min_i = -np.where(np.flip(np.array(L_s)/np.array(r_s)) >= L_r_min*0.99)[0][0]
		r_s, L_s = r_s[L_r_max_i:L_r_min_i], L_s[L_r_max_i:L_r_min_i]
		m_s = [self.m_Lr(L_s[i], r_s[i]) for i in range(len(r_s))]
		t_s = [self.t_r(r) for r in r_s]
		self.r, self.L, self.t, self.m = r_s, L_s, t_s, m_s
		return r_s, L_s, t_s, m_s

	def plot(self):
		#AT.plot_dual(self.r, self.L, self.m, "Radius [cm]", "Length [cm]", "Mass [kg]")
		AT.plot_single(np.array(self.L)/np.array(self.r), self.m, "Length/Radius", "Mass [kg]")

	def select(self, m=None, L_r=None, SF=1.0):
		if m is not None:
			selec = np.where(np.array(self.m) >= m)[0][0]
		elif L_r is not None:
			selec = np.where(np.array(self.L)/np.array(self.r) <= L_r)[0][0]
		else:
			raise ValueError("Either m or L_r should be specified")
		return self.r[selec], self.L[selec], self.t[selec]*SF, self.m[selec]*SF

	def V(self, L, r):
		return np.pi * r ** 2 * L + 4 / 3 * np.pi * r ** 3

	def L_r(self, r):
		return -4 * r / 3 + self.V_req / (r ** 2 * np.pi)

	def m_Lr(self, L, r, t=None):
		if t is None:
			t = self.t_r(r)
		return (2 * np.pi * r * L * t + 4 * np.pi * r ** 2 * t) * self.rho

	def t_r(self, r):
		return r / self.s_max * self.P

#O2_tank = tank(0.01, 0.5, 0.003, 13e6)
#O2_tank.plot()
#r, L, t, m = O2_tank.select(m=0.825)
#print(r, L, t, m)

#He_tank = tank(0.01, 1, 0.421, 30e6)
#He_tank.plot()
#r, L, t, m = He_tank.select(L_r=4.25)
#print(r, L, t, m)

#Hydrazine_tank = tank(0.01, 1, 0.82974, 8.75e6)
#Hydrazine_tank.plot()
#r, L, t, m = Hydrazine_tank.select(L_r=4.25)
#print(r, L, t, m)

#H2O2_tank = tank(0.01, 1, 1.0793, 8.75e6)
#H2O2_tank.plot()
#r, L, t, m = H2O2_tank.select(L_r=4.25)
#print(r, L, t, m)

#RCS_tank = tank(0.01, 1, 0.035, 8.75e6)
#RCS_tank.plot()
#r, L, t, m = RCS_tank.select(L_r=4)
#print(r, L, t, m)

RCS_tank = tank(0.01, 1, 0.21, 8.75e6)
#RCS_tank.plot()
r, L, t, m = RCS_tank.select(L_r=4)
print(r, L, t, m)