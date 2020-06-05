import numpy as np
import sympy as smp
from scipy.optimize import least_squares
import copy
from matplotlib import pyplot as plt
import EKF_class as EKF
from mpl_toolkits.mplot3d import Axes3D
from csv_impo import dataset
import math

p_l = [np.array([0.500,25,25]),
       np.array([0.1,-25,25]),
       np.array([0.1,25,-25]),
       np.array([0,-25,-25])]

r_l = smp.symbols('r0:%d'%len(p_l))
e_l = smp.symbols('e0:%d'%len(p_l))
b_l = smp.symbols('b0:%d'%len(p_l))
v_l = smp.symbols('v0:%d'%len(p_l))

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

xv= (x,y,z)
vv= (u,v,w)
av= (ax,ay,az)
eval_list = xv+vv+av+tuple([t])
data = dataset()
dt = 1
T = data.t[-1]
x0 = np.concatenate([data.r_p[0], data.v0, data.a0,np.array([0])])

sig = [3/20*10**-3 for r in r_l]+[0.05*np.pi/180 for e in e_l]+[0.05*np.pi/180  for b in b_l]+[1/90*10**-3 for v in v_l]+[(3.5e-6)*9.81 for a in av]+[1e-9]
sigq = [dt for r in r_l]+[0.001 for e in e_l]+[0.001 for b in b_l]+[0.0001 for v in v_l]+[0 for a in av]+[0]
h_val = 1e-3

dim = len(eval_list)
r_l += e_l + b_l + v_l + av + tuple([t])

def sign(input):
    if x<0:
        return -1
    else :
        return 1



p_l = [np.dot(data.trans_mat,(p + data.R0)) for p in p_l]

def af(i,time):
    dat = data.acc[i]
    i_p = math.floor(time/data.dt)
    i_n = math.ceil(time/data.dt)
    return (dat[i_n]-dat[i_p])/data.dt*(time-data.t[i_p])+dat[i_p]

r = EKF.vector_func([smp.sqrt(np.dot((xv-p),(xv-p))) for p in p_l] + [smp.asin((xv-p)[2]/smp.sqrt(np.dot((xv-p),(xv-p)))) for p in p_l] + [smp.acos((xv-p)[0]/smp.sqrt(np.dot((xv-p),(xv-p))))*smp.sign((xv-p)[1]) for p in p_l]  + [np.dot((vv),(xv-p))/smp.sqrt(np.dot((xv-p),(xv-p))) for p in p_l] + [a for a in av]+[t],eval_list)
r.lamb()
print(r.eval(x0))
xi0 = EKF.vector_func([x+u*dt+ax*dt**2/2,y+v*dt+ay*dt**2/2,z+w*dt+az*dt**2/2,u+ax*dt,v+ay*dt,w+az*dt,0,0,0,t+dt],eval_list) #implement trjaectory here
xi0.lamb()
xi0.main[6] = lambda x,y,z,u,v,w,ax,ay,az,t_temp : af(0,t_temp)
xi0.main[7] = lambda x,y,z,u,v,w,ax,ay,az,t_temp : af(1,t_temp)
xi0.main[8] = lambda x,y,z,u,v,w,ax,ay,az,t_temp : af(2,t_temp)



def ri1(r_val,r_per):
    def r_min(x):
        res_vect = r.eval(x)-r_val
        r_v = res_vect[:len(p_l)]
        e_v = np.multiply(res_vect[len(p_l):len(p_l)*2],r_v)
        b_v = np.multiply(res_vect[len(p_l)*2:len(p_l)*3],r_v)
        return np.concatenate([r_v,e_v,b_v])
    res = least_squares(r_min,r_per)
    return r.eval(xi0.eval(res.x))

h = EKF.vector_func(r_l,r_l)
Q = np.diag(sigq)*0
R = np.diag(sig)
r0 = r.eval(x0)


filter = EKF.EKF(ri1,h,Q,R,r0,comp = "num",h_val=h_val,n=len(r_l))

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
        z_func = r.eval(x_ref) + np.random.normal(0,sig,len(r_l))
        filter.corr(z_func,x_per)
    
    #Least square position
    def r_min(x):
            return r.eval(x)-filter.x

    res = least_squares(r_min,x_per)
    x_per = res.x
    l_ref[n] = x_ref
    l_out[n] = x_per


r = data.R0[0]
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0.0:pi/2:20j, 0.0:pi/2:20j]
xp = r*sin(phi)*cos(theta)
yp = r*sin(phi)*sin(theta)
zp = r*cos(phi)


plt.figure("trajectory")
ax = plt.axes(projection='3d')
ax.scatter3D(np.transpose(l_ref)[0],np.transpose(l_ref)[1],np.transpose(l_ref)[2])
ax.scatter3D(np.transpose(data.r_p)[0],np.transpose(data.r_p)[1],np.transpose(data.r_p)[2])
ax.scatter3D(np.transpose(l_out)[0],np.transpose(l_out)[1],np.transpose(l_out)[2])
ax.plot_surface(xp, yp, zp,  rstride=1, cstride=1, color='c', alpha=0.6, linewidth=0)
ax.scatter3D(np.transpose(p_l)[0],np.transpose(p_l)[1],np.transpose(p_l)[2])
ax.scatter3D(data.mars_r[0],data.mars_r[1],data.mars_r[0])
ax.set_xlim([0,4000])
ax.set_ylim([0,4000])
ax.set_zlim([0,4000])


plt.figure("deviation")
ax = plt.axes()
ax.plot(np.abs(np.transpose(l_ref)[0]-np.transpose(l_out)[0]))
ax.plot(np.abs(np.transpose(l_ref)[1]-np.transpose(l_out)[1]))
ax.plot(np.abs(np.transpose(l_ref)[2]-np.transpose(l_out)[2]))
plt.show()
