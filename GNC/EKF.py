import numpy as np
import sympy as smp
from sympy.solvers import solve

x = smp.Symbol('x')
y = smp.Symbol('y')

r1 = smp.Symbol('r1')
r2 = smp.Symbol('r2')
r3 = smp.Symbol('r3')

r_l = np.array([r1,r2,r3])

xb = np.array([x,y])
p_l = [np.array([3,1]),
    np.array([1,3]),
    np.array([1,1])]

r_sys = np.array([-r**2 for r in r_l])
r_sys += np.array([np.dot((xb-p),(xb-p)) for p in p_l])

u = np.array(solve(r_sys,xb.tolist()+[r1])[:-1])


