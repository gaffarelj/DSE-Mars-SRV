import numpy as np
import sympy as smp
from matplotlib import pyplot as plt
import EKF_class as EKF
from mpl_toolkits.mplot3d import Axes3D
import math
from csv_orbit import *
def dev(plot=false):
    data = dataset()

    dt = 0.1
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

    x0 = np.concatenate([data.r_p[0],data.v0,data.a0,[0]])
    h = EKF.vector_func(x_l,x_l)
    Q = np.diag(sig_m)*0
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
    def xf(x_in,xig):
        return xi0.eval(x_in)
    filter = EKF.EKF(xf,h,Q,R,x0,fcomp="num",h_val=1e-3,n=len(x_l))

    l_ref = np.ndarray([int(T/dt),dim])
    l_out = np.ndarray([int(T/dt),dim])

    x_per = x0
    x_ref = x0
    #print(axf(*np.concatenate([x0,np.array([0])])))
    for n in range(int(T/dt)):
        if plot:
            print(n)
            print(np.abs(x_per-x_ref))
        #update thuth
        x_ref = xi0.eval(x_ref)
        #forecast step
        filter.step(x_per)

        #correction step
        if True: 
            z_func = xi0.eval(x_ref) + np.random.normal(0,sig_a,len(x_l))
            filter.corr(z_func,x_per)
        
        x_per = filter.x
        l_ref[n] = x_ref
        l_out[n] = x_per

    if plot:
        print(filter.P)
        r = data.R0[0]
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:20j, 0.0:2*pi:20j]
        xp = r*sin(phi)*cos(theta)
        yp = r*sin(phi)*sin(theta)
        zp = r*cos(phi)

        plt.figure("trajectory")
        ax = plt.axes(projection='3d')
        ax.scatter3D(np.transpose(l_ref)[0],np.transpose(l_ref)[1],np.transpose(l_ref)[2])
        ax.scatter3D(np.transpose(data.r_p)[0],np.transpose(data.r_p)[1],np.transpose(data.r_p)[2])
        ax.scatter3D(np.transpose(l_out)[0],np.transpose(l_out)[1],np.transpose(l_out)[2])
        ax.plot_surface(xp, yp, zp,  rstride=1, cstride=1, color='c', alpha=0.6, linewidth=0)
        ax.scatter3D(data.mars_r[0],data.mars_r[1],data.mars_r[0])
        #ax.set_xlim([0,4000])
        #ax.set_ylim([0,4000])
        #ax.set_zlim([0,4000])


        plt.figure("deviation")
        ax = plt.subplot(121)
        ax.plot(np.arange(0,T-dt,dt),np.abs(np.transpose(l_ref)[0]-np.transpose(l_out)[0])*10**3,label="Downrange position deviation")
        ax.plot(np.arange(0,T-dt,dt),np.abs(np.transpose(l_ref)[1]-np.transpose(l_out)[1])*10**3,label="Crossrange position deviation")
        ax.plot(np.arange(0,T-dt,dt),np.abs(np.transpose(l_ref)[2]-np.transpose(l_out)[2])*10**3,label="Altitude position deviation")
        ax.legend()
        ax.set_ylabel("deviation position[m]")
        ax.set_xlabel("aproach time [s]")
        ax = plt.subplot(122)
        ax.plot(np.arange(0,T-dt,dt),np.abs(np.transpose(l_ref)[3]-np.transpose(l_out)[3])*10**3,label="Downrange velocity deviation")
        ax.plot(np.arange(0,T-dt,dt),np.abs(np.transpose(l_ref)[4]-np.transpose(l_out)[4])*10**3,label="Crossrange velocity deviation")
        ax.plot(np.arange(0,T-dt,dt),np.abs(np.transpose(l_ref)[5]-np.transpose(l_out)[5])*10**3,label="Altitude velocity deviation")
        ax.legend()
        ax.set_ylabel("deviation velocity [m/s]")
        ax.set_xlabel("aproach time [s]")
        plt.show()

return (l_ref[-1]-l_out[-1])[:6]*10**3