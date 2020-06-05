import numpy as np
import sympy as smp
from matplotlib import pyplot as plt
import EKF_class as EKF
from mpl_toolkits.mplot3d import Axes3D
import math


data = dataset()

dt = 1
T = data.t[-1]

x = smp.Symbol('x')
y = smp.Symbol('y')
z = smp.Symbol('z')
u = smp.Symbol('u')
v = smp.Symbol('v')
w = smp.Symbol('w')
ax = smp.Symbol('ax')
ay = smp.Symbol('ay')
az = smp.Symbol('az')

t = smp.Symbol('t')

p = (x,y,z)
pd1 = (u,v,w)
pd2 = (ax,ay,az)

x_l = p + pd1 + pd2 + tuple([t])

dim = len(x_l)

sig_m = np.array([ 1 for var in p] + [1 for var in pd1] + [(3.5e-6)*9.81 for var in pd2] + [1e-9])
sig_a = np.array([ 0 for var in p] + [0 for var in pd1] + [(3.5e-6)*9.81 for var in pd2] + [1e-9])

x0 = 
h = EKF.vector_func(x_l,x_l)
Q = np.diag(sig)*0
R = np.diag(sig_m)

def af(i,time):
    dat = data.acc[i]
    i_p = math.floor(time/data.dt)
    i_n = math.ceil(time/data.dt)
    return (dat[i_n]-dat[i_p])/data.dt*(time-data.t[i_p])+dat[i_p]

xi0 = EKF.vector_func([x+u*dt+ax*dt**2/2,y+v*dt+ay*dt**2/2,z+w*dt+az*dt**2/2,u+ax*dt,v+ay*dt,w+az*dt,0,0,0,t+dt],x_l) 
xi0.lamb()
xi0.main[6] = lambda x,y,z,u,v,w,ax,ay,az,t_temp : af(0,t_temp)
xi0.main[7] = lambda x,y,z,u,v,w,ax,ay,az,t_temp : af(1,t_temp)
xi0.main[8] = lambda x,y,z,u,v,w,ax,ay,az,t_temp : af(2,t_temp)

filter = EKF.EKF(xi0,h,Q,R,x0,comp = "num",h_val=h_val,n=len(x_l))

l_ref = np.ndarray([int(T/dt),dim])
l_out = np.ndarray([int(T/dt),dim])

x_per = x0
x_ref = x0
#print(axf(*np.concatenate([x0,np.array([0])])))
for n in range(int(T/dt)):
    print(n)
    print(np.abs(x_per-x_ref))
    #update thuth
    x_ref = xi0.eval(x_ref)
    #forecast step
    filter.step(x_per)

    #correction step
    if True: 
        z_func = r.eval(x_ref) + np.random.normal(0,sig_a,len(r_l))
        filter.corr(z_func,x_per)
    
    x_per = filter.x
    l_ref[n] = x_ref
    l_out[n] = x_per

