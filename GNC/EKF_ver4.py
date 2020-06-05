import numpy as np
import sympy as smp
from scipy.optimize import least_squares
import copy
from matplotlib import pyplot as plt
import EKF_class as EKF
from mpl_toolkits.mplot3d import Axes3D

p_l = [np.array([-10,-5,0.500]),
       np.array([100,-20,0.020]),
       np.array([100,20,0.300]),
       np.array([-10,5,1.200])]
    
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

xv= (x,y,z)
vv= (u,v,w)
av= (ax,ay,az)

dt = 0.1
T = 100
x0 = np.array([500,0,80,-5,0,-0.800,0,0,0])

sig = [3/20*10**-3 for r in r_l]+[0.05*np.pi/180 for e in e_l]+[0.05*np.pi/180  for b in b_l]+[1/90*10**-3 for v in v_l]+[(3.5e-6)*9.81 for a in av]
h_val = 1e-3

dim = len(xv+vv+av)
r_l += e_l + b_l + v_l + av

def sign(input):
    if x<0:
        return -1
    else :
        return 1

r = EKF.vector_func([smp.sqrt(np.dot((xv-p),(xv-p))) for p in p_l] + [smp.asin((xv-p)[2]/smp.sqrt(np.dot((xv-p),(xv-p)))) for p in p_l] + [smp.acos((xv-p)[0]/smp.sqrt(np.dot((xv-p),(xv-p))))*smp.sign((xv-p)[1]) for p in p_l]  + [np.dot((vv),(xv-p))/smp.sqrt(np.dot((xv-p),(xv-p))) for p in p_l] + [a for a in av],xv+vv+av)
r.lamb()
print(r.eval(x0))
xi0 = EKF.vector_func([x+u*dt+ax*dt**2/2,y+v*dt+ay*dt**2/2,z+w*dt+az*dt**2/2,u+ax*dt,v+ay*dt,w+az*dt,ax,ay,az],xv+vv+av) #implement trjaectory here
xi0.lamb()

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
Q = np.eye(len(r_l))*0
R = np.diag(sig)
r0 = r.eval(x0)


filter = EKF.EKF(ri1,h,Q,R,r0,comp = "num",h_val=h_val,n=len(r_l))

l_ref = np.ndarray([int(T/dt),dim])
l_out = np.ndarray([int(T/dt),dim])

x_per = x0
x_ref = x0

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


plt.figure("x,y")
ax = plt.axes(projection='3d')
ax.scatter3D(np.transpose(l_ref)[0],np.transpose(l_ref)[1],np.transpose(l_ref)[2])
ax.scatter3D(np.transpose(l_out)[0],np.transpose(l_out)[1],np.transpose(l_out)[2])
ax.scatter3D(np.transpose(p_l)[0],np.transpose(p_l)[1],np.transpose(p_l)[2])
ax.scatter3D([0],[0],[0])
plt.figure("deviation")
ax = plt.axes()
ax.plot(np.abs(np.transpose(l_ref)[0]-np.transpose(l_out)[0]))
ax.plot(np.abs(np.transpose(l_ref)[1]-np.transpose(l_out)[1]))
ax.plot(np.abs(np.transpose(l_ref)[2]-np.transpose(l_out)[2]))
plt.show()
