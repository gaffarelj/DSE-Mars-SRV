import numpy as np
import sympy as smp
from scipy.optimize import least_squares
import copy
from matplotlib import pyplot as plt
import EKF_class as EKF


r1 = smp.Symbol('r1')
r2 = smp.Symbol('r2')
r3 = smp.Symbol('r3')
v1 = smp.Symbol('v1')
v2 = smp.Symbol('v2')
v3 = smp.Symbol('v3')

x = smp.Symbol('x')
y = smp.Symbol('y')
u = smp.Symbol('u')
v = smp.Symbol('v')

xv= [x,y]
vv= [u,v]


r_l = [r1,r2,r3,v1,v2,v3]



p_l = [np.array([3,1]),
       np.array([1,3]),
       np.array([1,1])]



#x of shape [x,y,u,v]

dt = 0.01
T = 100
x0 = np.array([2,2,0.1,0.1])
dim = 4
sig = 0.005



r = EKF.vector_func([smp.sqrt(np.dot((xv-p),(xv-p))) for p in p_l] + [np.dot((vv),(xv-p))/smp.sqrt(np.dot((xv-p),(xv-p))) for p in p_l],xv+vv)
r.lamb()
"""
xi1 = EKF.vector_func([-r1**2/18+r3**2/18+11/2,
                      (-r2**2+r3**2+8)/4,
                      9*(v2*r2*(r2**2-r3**2+4)+v3*r3*(r3**2-r2**2-12)/(4*(r1**2-r3**2-81))),
                      -v2*r2/2+v3*r3/2],
                      r_l)
"""

xi1 = EKF.vector_func([(-r1**2+r3**2+8)/4,
                      (-r2**2+r3**2+8)/4,
                      -v1*r1/2+v3*r3/2,
                      -v2*r2/2+v3*r3/2],
                      r_l)


xi0 = EKF.vector_func([x+u*dt,y+v*dt,u,v],xv+vv) #implement trjaectory here
xi0r = xi0.subs_func(xi1)
xi0.lamb()
xi0r.var = r_l
xi0r.lamb()
r_sys = r.subs_func(xi0r,xv+vv)
r_sys.var = r_l

h = EKF.vector_func(r_l,r_l)


Q = np.eye(6)*0
R = np.eye(6)*sig
r0 = r.eval(x0)


filter = EKF.EKF(r_sys,h,Q,R,r0)

l_ref = np.ndarray([int(T/dt),dim])
l_out = np.ndarray([int(T/dt),dim])

x_per = x0
x_ref = x0
print(xi0.eval(x0))
print(r_sys.eval(r0))

for n in range(int(T/dt)):
       #print(n)
       #update thuth
       x_ref = xi0.eval(x_ref)
       #forecast step
       filter.step()
       
       #correction step
       if n%1==0:
              z = r.eval(x_ref) + np.random.normal(0,sig,len(r_l))       
              filter.corr(z)
       #print(filter.x[0]**2-filter.x[2]**2-81)
       #Least square position
       def r_min(x):
              return np.array(r.eval(x))-np.array(filter.x)
       
       res = least_squares(r_min,x_per)
       x_per = res.x

       l_ref[n] = x_ref
       l_out[n] = x_per


plt.figure("x,y")
plt.scatter(np.transpose(l_ref)[0],np.transpose(l_ref)[1])
plt.scatter(np.transpose(l_out)[0],np.transpose(l_out)[1])
plt.scatter(np.transpose(p_l)[0],np.transpose(p_l)[1])
plt.figure("deviation")
N = 200
plt.plot(np.arange(0,T,dt),np.convolve(np.transpose(np.abs(l_ref-l_out))[0], np.ones((N,))/N, mode='same'))
plt.plot(np.arange(0,T,dt),np.convolve(np.transpose(np.abs(l_ref-l_out))[1], np.ones((N,))/N, mode='same'))
plt.show()


#filter.corr(np.array([2,2,3]))
