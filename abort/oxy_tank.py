import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT

def V(L, r):
	return np.pi * r ** 2 * L + 4 / 3 * np.pi * r ** 3

def L_r(r):
	return -4 * r / 3 + 3 / (1000 * r ** 2 * np.pi)

def m(L, r, t=None, rho=2700):
	if t is None:
		t = t_r(r)
	return (2 * np.pi * r * L * t + 4 * np.pi * r ** 2 * t) * rho

def t_r(r):
	return r / 21.23

r_s = np.arange(0.01, 0.5, 0.001)
L_s = [L_r(r) for r in r_s]
m_s = [m(L_s[i], r_s[i]) for i in range(len(r_s))]

L_neg = np.where(np.array(L_s) <= 0)[0][0]

AT.plot_dual([r*100 for r in r_s[:L_neg]], [L*100 for L in L_s[:L_neg]], m_s[:L_neg], "Radius [cm]", "Length [cm]", "Mass [kg]")

m_80 = np.where(np.array(m_s) >= 0.8)[0][0]
print(m_s[m_80], L_s[m_80], r_s[m_80], t_r(r_s[m_80]))