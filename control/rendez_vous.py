import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act

import numpy as np
import sympy as sp
from matplotlib import pyplot as plt


#=====================================================================================================================================================================================================================
# Function to compute Thrust 
#=====================================================================================================================================================================================================================
def vac_thrust(DeltaV,Isp,Mbegin,tb,De=0,pe=0):
    """computes vacuum thrust from DeltaV. 
        Mbegin=mass at the start of the maneuver, tb=burn time, De=exit diameter of engine/thruster, pe=exhaust exit pressure of engine/thruster
    """
    Ae=np.pi/4*De*De
    thrust=(np.exp(DeltaV/(Isp*9.80665))*Mbegin-Mbegin)/(tb*np.exp(DeltaV/(Isp*9.80665)))*Isp*9.80665+Ae*(pe)
    return thrust

#=====================================================================================================================================================================================================================
#Node properties
#=====================================================================================================================================================================================================================
mu=0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
R=3389.5*10**3                  #[m] volumetric mean radius
h_node=500*10**3                #[m]

period = np.sqrt((R+h_node) ** 3 / mu)
omega  = period / (2 * np.pi)

#=====================================================================================================================================================================================================================
#Vehicle properties
#=====================================================================================================================================================================================================================

#Implement Mass of Charon at t0!
m0   =53248.16053543461         #previous value was: 100000
Isp  = 221                      #previous value was: 390

#=====================================================================================================================================================================================================================
#Vbar approach
#=====================================================================================================================================================================================================================
#Proximity operations A:
x0_A = -1000         #m
x1_A = -250          #m
t_A  = 5 * period   #s
Vx_A   = (x0_A - x1_A) / t_A #m/s
deltaV_A_0 = deltaV_A_1 = Vx_A

#Proximity operations B:
x0_B = x1_A         #m
x1_B = -30           #m
t_B  = 60 * 60      #s
Vx_B  = (x0_B - x1_B) / t_B #m/s
deltaV_B_0 = deltaV_B_1 = Vx_B
#Docking:
x0_d = x1_B         #m
x1_d = -3            #m
t_d  = (5*((x0_d - x1_d)/10)) * 60       #m
Vx_d  = (x0_d - x1_d) / t_d #m/s
deltaV_d_0 = deltaV_d_1 = Vx_d

y0 = z0 = 0

deltaV_tot = deltaV_A_0 + deltaV_A_1 + deltaV_B_0 + deltaV_B_1 + deltaV_d_0 + deltaV_d_1

print('Velocities:',Vx_A,Vx_B,Vx_d)


#=====================================================================================================================================================================================================================
#Rbar approach
#=====================================================================================================================================================================================================================
# #Proximity operations A:
# z0_A = 1000         #m
# z1_A = 250          #m
# t_A  = 5 * period   #s
#
# #Proximity operations B:
# z0_B = x1_A         #m
# z1_B = 30           #m
# t_B  = 70 * 60      #s
# #Docking:
# z0_d = x1_B         #m
# z1_d = 3            #m
# t_d  = 5 * 60       #m
#
# x0 = y0 = 0

#=====================================================================================================================================================================================================================
#Navigation measurement errors
#=====================================================================================================================================================================================================================
error_xm = 1.
error_ym = 1.
error_zm = 1.

error_xdot = 0.1
error_ydot = 0.1
error_zdot = 0.1

def error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t):
    delta_x = error_xm + 6 * error_zm * (omega * t - np.sin(omega * t)) +  error_xdot * (4 / omega * np.sin(omega * t) - 3 * t) + 2 / omega * error_zdot * (1 - np.cos(omega * t))
    return delta_x
def error_y(error_ym,omega,t):
    delta_y = error_ym * np.cos(omega * t) + 1 / omega * error_ydot * np.sin (omega * t)
    return delta_y
def error_z(error_xm,error_zm,error_zdot,omega,t):
    delta_z = error_zm * (4 - 3 * np.cos(omega * t)) + 2 / omega * error_xdot * (np.cos(omega * t) -1) + 1 / omega * error_zdot * np.sin(omega * t)
    return delta_z


#=====================================================================================================================================================================================================================
# Switches
#=====================================================================================================================================================================================================================

plotting=True       #Do you wanna plot? no=False

#=====================================================================================================================================================================================================================
# Hill equations of motion
#=====================================================================================================================================================================================================================

def thrust_x(x,z,omega,t):
    t = sp.Symbol('t')
    xdot = sp.diff(x)
    xdotdot = sp.diff(xdot)
    zdot = sp.diff(z)
    zdotdot = sp.diff(zdot)
    gamma_x = (xdotdot - 2 * omega * zdot)
    return gamma_x

def thrust_y(y,omega,t):
    t = sp.Symbol('t')
    ydot = sp.diff(y)
    ydotdot = sp.diff(ydot)
    gamma_y = (ydotdot + omega ** 2 * y)
    return gamma_y

def thrust_z(x,z,omega,t):
    t = sp.Symbol('t')
    xdot = sp.diff(x)
    zdot = sp.diff(z)
    zdotdot = sp.diff(zdot)
    gamma_z = (zdotdot + 2 * omega * xdot - 3 * omega ** 2 * z)
    return gamma_z

#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
dt = 0.01
t = 0.
mp = 0.
m = m0
delta_x = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
delta_y = error_y(error_ym,omega,t)
delta_z = error_z(error_xm,error_zm,error_zdot,omega,t)
x = x0_A + delta_x
y = y0 + delta_y
z = z0 + delta_z

fx0 = thrust_x(x,z,omega,t)
fy0 = thrust_y(y,omega,t)
fz0 = thrust_z(x,z,omega,t)
f_array   = np.array([[fx0,fy0,fz0]])
X_array   = np.array([[x,y,z]])
t_array   = np.array(t)
mp_array  = np.array(mp)

#Proximity operations A:
while x >= x0_A and x < x1_A:
    delta_x_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y_new = error_y(error_ym,omega,t)
    delta_z_new = error_z(error_xm,error_zm,error_zdot,omega,t)

    delta_x_rel = delta_x_new - delta_x
    delta_y_rel = delta_y_new - delta_y
    delta_z_rel = delta_z_new - delta_z
    #print('deltas: ', delta_x_rel,delta_y_rel,delta_z_rel)
    Vx = Vx_A
    x  += Vx*dt + delta_x_rel
    y  = delta_y_rel
    z  = delta_z_rel
    fx = m * thrust_x(x,z,omega,t)
    fy = m * thrust_y(y,omega,t)
    fz = m * thrust_z(x,z,omega,t)
    ftot = fx + fy + fz
    print("this is thrust",ftot)
    mp = act.RCSpropellant(ftot,dt,Isp)
    m -= mp
    t += dt


    f_array  = np.append(f_array,[[fx,fy,fz]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp)
    t_array  = np.append(t_array,t)

    delta_x = delta_x_new
    delta_y = delta_y_new
    delta_z = delta_z_new

    print('phase A')
    print('x: ',x)
    print(delta_x,delta_y,delta_z)
    print('============================')

while x >= x0_B and x < x1_B:
    delta_x_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y_new = error_y(error_ym,omega,t)
    delta_z_new = error_z(error_xm,error_zm,error_zdot,omega,t)

    delta_x_rel = delta_x_new - delta_x
    delta_y_rel = delta_y_new - delta_y
    delta_z_rel = delta_z_new - delta_z

    Vx = Vx_B
    x  += Vx*dt + delta_x_rel
    y  = delta_y_rel
    z  = delta_z_rel
    fx = m * thrust_x(x,z,omega,t)
    fy = m * thrust_y(y,omega,t)
    fz = m * thrust_z(x,z,omega,t)
    ftot = fx + fy + fz
    mp = act.RCSpropellant(ftot,dt,Isp)
    m -= mp
    t += dt

    # f_array = np.append(f_array,[[fx,fy,fz]],axis=0)
    # X_array = np.append(X_array,[[x,y,z]],axis=0)
    # mp_array = np.append(mp_array,mp)
    # t_array = np.append(t_array,t)

    delta_x = delta_x_new
    delta_y = delta_y_new
    delta_z = delta_z_new

    print('phase B')
    print('x: ',x)
    print('============================')

while x >= x0_d and x < x1_d:
    print('phase d')
    delta_x_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y_new = error_y(error_ym,omega,t)
    delta_z_new = error_z(error_xm,error_zm,error_zdot,omega,t)

    delta_x_rel = delta_x_new - delta_x
    delta_y_rel = delta_y_new - delta_y
    delta_z_rel = delta_z_new - delta_z
    Vx += Vx_B
    x  = Vx*dt + delta_x_rel
    y  = delta_y_rel
    z  = delta_z_rel
    fx = m * thrust_x(x,z,omega,t)
    fy = m * thrust_y(y,omega,t)
    fz = m * thrust_z(x,z,omega,t)
    ftot = fx + fy + fz
    mp = act.RCSpropellant(ftot,dt,Isp)
    m -= mp
    t += dt

    # f_array = np.append(f_array,[[fx,fy,fz]],axis=0)
    # X_array = np.append(X_array,[[x,y,z]],axis=0)
    # mp_array = np.append(mp_array,mp)
    # t_array = np.append(t_array,t)

    delta_x = delta_x_new
    delta_y = delta_y_new
    delta_z = delta_z_new

    print('phase d')
    print('x: ',x)
    print('============================')


print(f_array[:,1])

#=====================================================================================================================================================================================================================
# Thrust to initialize and/or stop a given phase.
# ASSUMPTION: thrusters fire for 0.5 sec
#=====================================================================================================================================================================================================================



#=====================================================================================================================================================================================================================
# Plotting
#=====================================================================================================================================================================================================================
if plotting:
    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    fig.suptitle('Relative position wrt target',fontsize=16)
    #Fx vs time
    axs[0][0].plot(t_array,X_array[:,0],color="navy")
    axs[0][0].grid(color="gainsboro")
    axs[0][0].set_xlabel("Time [s]")
    axs[0][0].set_ylabel("X position [m]")

    #Fx vs time
    axs[0][1].plot(t_array,X_array[:,1],color="navy")
    axs[0][1].grid(color="gainsboro")
    axs[0][1].set_xlabel("Time [s]")
    axs[0][1].set_ylabel("Y position [m]")

    #Fx vs time
    axs[0][2].plot(t_array,X_array[:,2],color="navy")
    axs[0][2].grid(color="gainsboro")
    axs[0][2].set_xlabel("Time [s]")
    axs[0][2].set_ylabel("Z postion [m]")
    plt.show()

    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    fig.suptitle('RCS thrust',fontsize=16)
    #Fx vs time
    axs[0][0].plot(t_array,f_array[:,0],color="navy")
    axs[0][0].grid(color="gainsboro")
    axs[0][0].set_xlabel("Time [s]")
    axs[0][0].set_ylabel("Thrust [N]")

    #Fx vs time
    axs[0][1].plot(t_array,f_array[:,1],color="navy")
    axs[0][1].grid(color="gainsboro")
    axs[0][1].set_xlabel("Time [s]")
    axs[0][1].set_ylabel("Thrust [N]")

    #Fx vs time
    axs[0][2].plot(t_array,f_array[:,2],color="navy")
    axs[0][2].grid(color="gainsboro")
    axs[0][2].set_xlabel("Time [s]")
    axs[0][2].set_ylabel("Thrust [N]")
    plt.show()
