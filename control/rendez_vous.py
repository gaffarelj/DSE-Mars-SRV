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
    Mprop=(np.exp(DeltaV/(Isp*9.80665))*Mbegin-Mbegin)/(np.exp(DeltaV/(Isp*9.80665)))
    return thrust, Mprop

#=====================================================================================================================================================================================================================
#Node properties
#=====================================================================================================================================================================================================================
<<<<<<< HEAD
mu=0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
R=3389.5*10**3                  #[m] volumetric mean radius
h_node=500*10**3                #[m]
h_phasing=609.74*10**3
period =2*np.pi* np.sqrt((R+h_node) ** 3 / mu)
time_hohmann=np.pi*np.sqrt(((h_node+h_phasing+2*R)/2)**3*1/mu)
=======
mu     = 0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
R      = 3389.5*10**3                  #[m] volumetric mean radius
h_node = 500*10**3
g_earth= 9.81             #[m]

period =  2*np.pi* np.sqrt((R+h_node) ** 3 / mu)
>>>>>>> 96efb5901be3b77f1984cbaf043bf3b6990e77a8
omega  = (2 * np.pi)/period

#=====================================================================================================================================================================================================================
#Vehicle properties
#=====================================================================================================================================================================================================================

#Implement Mass of Charon at t0!
m0   =53248.16053543461         #previous value was: 100000
Isp  = 221                      #previous value was: 390
Ix = 2875350.278 #kg/m^2
Iy = 306700.3372 #kg/m^2
Iz = 2875350.278 #kg/m^2

#=====================================================================================================================================================================================================================
#Vbar approach properties
#=====================================================================================================================================================================================================================
#Proximity operations A:
x0_A = -1000         #m
x1_A = -250          #m
t_A  =  period   #s
Vx_A   = (x1_A - x0_A) / t_A #m/s
deltaV_A_0 = deltaV_A_1 = Vx_A

#Proximity operations B:
x0_B = x1_A         #m
x1_B = -30           #m
t_B  = 60 * 60      #s
Vx_B  = (x1_B - x0_B) / t_B #m/s
deltaV_B_0 = deltaV_B_1 = Vx_B
#Docking:
x0_d = x1_B         #m
x1_d = -3            #m
t_d  = (5*((x1_d - x0_d)/10)) * 60       #m
Vx_d  = (x1_d - x0_d) / t_d #m/s
deltaV_d_0 = deltaV_d_1 = Vx_d

y0 = z0 = 0

deltaV_tot = deltaV_A_0 + deltaV_A_1 + deltaV_B_0 + deltaV_B_1 + deltaV_d_0 + deltaV_d_1


#=====================================================================================================================================================================================================================
#Navigation measurement errors
#=====================================================================================================================================================================================================================
error_xm = 0.1
error_ym = 0.1
error_zm = 0.1

error_xdot = 0.01
error_ydot = 0.01
error_zdot = 0.01

def error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t):
    delta_x = error_xm + 6 * error_zm * (omega * t - np.sin(omega * t)) +  error_xdot * (4 / omega * np.sin(omega * t) - 3 * t) + 2 / omega * error_zdot * (1 - np.cos(omega * t))
    delta_xdot = 6 * error_zm * (omega - omega * np.cos(omega * t)) +  error_xdot * (4 * np.cos(omega * t) - 3) + 2 / omega * error_zdot * (omega *  np.sin(omega * t))
    delta_xdotdot = 6 * error_zm * (omega ** 2 * np.sin(omega * t)) +  error_xdot * (-4 * omega * np.sin(omega * t)) + 2 / omega * error_zdot * (omega ** 2 * np.cos(omega * t))

    return delta_x, delta_xdot, delta_xdotdot

def error_y(error_ym,omega,t):
    delta_y = error_ym * np.cos(omega * t) + 1 / omega * error_ydot * np.sin(omega * t)
    delta_ydot = -omega * error_ym * np.sin(omega * t) + error_ydot * np.cos(omega * t)
    delta_ydotdot = -omega ** 2 * error_ym * np.cos(omega * t) - omega * error_ydot * np.sin(omega * t)

    return delta_y, delta_ydot, delta_ydotdot

def error_z(error_xm,error_zm,error_zdot,omega,t):
    delta_z = error_zm * (4 - 3 * np.cos(omega * t)) + 2 / omega * error_xdot * (np.cos(omega * t) -1) + 1 / omega * error_zdot * np.sin(omega * t)
    delta_zdot=3*error_zm*omega*np.sin(omega*t)-2*error_xdot*np.sin(omega*t)+error_zdot*np.cos(omega*t)
    delta_zdotdot=3*error_zm*omega*omega*np.cos(omega*t)-2*error_xdot*omega*np.cos(omega*t)-error_zdot*omega*np.sin(omega*t)

    return delta_z, delta_zdot, delta_zdotdot


#=====================================================================================================================================================================================================================
# Switches
#=====================================================================================================================================================================================================================

plotting=True       #Do you wanna plot? no=False

#=====================================================================================================================================================================================================================
# Hill equations of motion
#=====================================================================================================================================================================================================================

def thrust_x(x,zdot,xdotdot,omega,t):
    gamma_x = (xdotdot - 2 * omega * zdot)
    return gamma_x

def thrust_y(y,ydotdot,omega,t):
    gamma_y = (ydotdot + omega ** 2 * y)
    return gamma_y

def thrust_z(z,zdotdot,xdot,omega,t):
    gamma_z = (zdotdot + 2 * omega * xdot - 3 * omega ** 2 * z)
    return gamma_z

#=====================================================================================================================================================================================================================
# Simulation
#=====================================================================================================================================================================================================================
dt = 1.
t  = 0.1
ta = tb = td = 0.1
tburn = 1.
mp = 0.
m = m0
m_deltaV = [m]

delta_x, delta_xdot, delta_xdotdot = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
delta_y, delta_ydot, delta_ydotdot = error_y(error_ym,omega,t)
delta_z, delta_zdot, delta_zdotdot = error_z(error_xm,error_zm,error_zdot,omega,t)
x   = x0_A + delta_x
xdot    = delta_xdot
xdotdot = delta_xdotdot
y   = y0 + delta_y
ydot    = delta_ydot
ydotdot = delta_ydotdot
z   = z0 + delta_z
zdot    = delta_zdot
zdotdot = delta_zdotdot
print('initial xyz: ', x,y,z)
print('initial delta xyz: ', delta_x,delta_y,delta_z)

fx0 = m * thrust_x(x,zdot,xdotdot,omega,t)
fy0 = m * thrust_y(y,ydotdot,omega,t)
fz0 = m * thrust_z(z,zdotdot,xdot,omega,t)

f_array   = np.array([[fx0,fy0,fz0]])
X_array   = np.array([[x,y,z]])
t_array   = np.array(t)
mp_array  = np.array(mp)

#DeltaV maneuver 1
t += tburn
thrust_deltaV1, mp_deltaV1 = vac_thrust(deltaV_A_0,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV1
f_array = np.append(f_array,[[thrust_deltaV1,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV1)
t_array  = np.append(t_array,t)

#Proximity operations A:
while x >= x0_A-1 and x < x1_A:
    Vx = Vx_A
    t += dt
    ta += dt
    delta_x_new, delta_xdot_new, delta_xdotdot_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y_new, delta_ydot_new, delta_ydotdot_new = error_y(error_ym,omega,t)
    delta_z_new, delta_zdot_new, delta_zdotdot_new = error_z(error_xm,error_zm,error_zdot,omega,t)

    delta_x_rel       = delta_x_new - delta_x
    delta_xdot_rel    = delta_xdot_new - delta_xdot
    delta_xdotdot_rel = delta_xdotdot_new - delta_xdotdot

    delta_y_rel       = delta_y_new - delta_y
    delta_ydot_rel    = delta_ydot_new - delta_ydot
    delta_ydotdot_rel = delta_ydotdot_new - delta_ydotdot

    delta_z_rel       = delta_z_new - delta_z
    delta_zdot_rel    = delta_zdot_new - delta_zdot
    delta_zdotdot_rel = delta_zdotdot_new - delta_zdotdot

    x       = x0_A + Vx*ta + delta_x_rel
    xdot    = Vx + delta_xdot_rel
    xdotdot = delta_xdotdot_rel

    y       = delta_y_rel
    ydot    = delta_ydot_rel
    ydotdot = delta_ydotdot_rel

    z       = delta_z_rel
    zdot    = delta_zdot_rel
    zdotdot = delta_zdotdot_rel

    fx = m * thrust_x(x,zdot,xdotdot,omega,t)
    fy = m * thrust_y(y,ydotdot,omega,t)
    fz = m * thrust_z(z,zdotdot,xdot,omega,t)
    ftot = abs(fx + fy + fz)
    mp = act.RCSpropellant(ftot,dt,Isp)
    m -= mp

    f_array  = np.append(f_array,[[fx,fy,fz]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp)
    t_array  = np.append(t_array,t)

    delta_x       = delta_x_new
    delta_xdot    = delta_xdot_new
    delta_xdotdot = delta_xdotdot_new

    delta_y       = delta_y_new
    delta_ydot    = delta_ydot_new
    delta_ydotdot = delta_ydotdot_new

    delta_z       = delta_z_new
    delta_zdot    = delta_zdot_new
    delta_zdotdot = delta_zdotdot_new

#Delta V maneuver 2
t += tburn
thrust_deltaV2, mp_deltaV2 = vac_thrust(deltaV_A_1,Isp,m,tburn,De=0,pe=0)
m -= mp_deltaV2
f_array = np.append(f_array,[[thrust_deltaV2,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV2)
t_array  = np.append(t_array,t)

#Delta V maneuver 3
t += tburn
thrust_deltaV3, mp_deltaV3 = vac_thrust(deltaV_B_0,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV3
f_array = np.append(f_array,[[thrust_deltaV3,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV3)
t_array  = np.append(t_array,t)

#Proximity operations B
while x >= x0_B-1 and x < x1_B:
    tb += dt
    Vx = Vx_B

    delta_x_new, delta_xdot_new, delta_xdotdot_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y_new, delta_ydot_new, delta_ydotdot_new = error_y(error_ym,omega,t)
    delta_z_new, delta_zdot_new, delta_zdotdot_new = error_z(error_xm,error_zm,error_zdot,omega,t)

    delta_x_rel       = delta_x_new - delta_x
    delta_xdot_rel    = delta_xdot_new - delta_xdot
    delta_xdotdot_rel = delta_xdotdot_new - delta_xdotdot

    delta_y_rel       = delta_y_new - delta_y
    delta_ydot_rel    = delta_ydot_new - delta_ydot
    delta_ydotdot_rel = delta_ydotdot_new - delta_ydotdot

    delta_z_rel       = delta_z_new - delta_z
    delta_zdot_rel    = delta_zdot_new - delta_zdot
    delta_zdotdot_rel = delta_zdotdot_new - delta_zdotdot


    x       = x0_B + Vx*tb + delta_x_rel
    xdot    = Vx + delta_xdot_rel
    xdotdot = delta_xdotdot_rel

    y       = delta_y_rel
    ydot    = delta_ydot_rel
    ydotdot = delta_ydotdot_rel

    z       = delta_z_rel
    zdot    = delta_zdot_rel
    zdotdot = delta_zdotdot_rel

    fx = m * thrust_x(x,zdot,xdotdot,omega,t)
    fy = m * thrust_y(y,ydotdot,omega,t)
    fz = m * thrust_z(z,zdotdot,xdot,omega,t)
    ftot = abs(fx + fy + fz)
    mp = act.RCSpropellant(ftot,dt,Isp)
    m -= mp


    f_array = np.append(f_array,[[fx,fy,fz]],axis=0)
    X_array = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp)
    t_array = np.append(t_array,t)

    delta_x       = delta_x_new
    delta_xdot    = delta_xdot_new
    delta_xdotdot = delta_xdotdot_new

    delta_y       = delta_y_new
    delta_ydot    = delta_ydot_new
    delta_ydotdot = delta_ydotdot_new

    delta_z       = delta_z_new
    delta_zdot    = delta_zdot_new
    delta_zdotdot = delta_zdotdot_new

    t  += dt

#Delta V maneuver 4
t += tburn
thrust_deltaV4, mp_deltaV4 = vac_thrust(deltaV_B_1,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV4
f_array = np.append(f_array,[[thrust_deltaV4,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV4)
t_array  = np.append(t_array,t)

#Delta V maneuver 5
t += tburn
thrust_deltaV5, mp_deltaV5 = vac_thrust(deltaV_d_0,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV5
f_array = np.append(f_array,[[thrust_deltaV5,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV5)
t_array  = np.append(t_array,t)

while x >= x0_d-1 and x < x1_d:
    Vx  = Vx_d
    td += dt

    delta_x_new, delta_xdot_new, delta_xdotdot_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y_new, delta_ydot_new, delta_ydotdot_new = error_y(error_ym,omega,t)
    delta_z_new, delta_zdot_new, delta_zdotdot_new = error_z(error_xm,error_zm,error_zdot,omega,t)

    delta_x_rel       = delta_x_new - delta_x
    delta_xdot_rel    = delta_xdot_new - delta_xdot
    delta_xdotdot_rel = delta_xdotdot_new - delta_xdotdot

    delta_y_rel       = delta_y_new - delta_y
    delta_ydot_rel    = delta_ydot_new - delta_ydot
    delta_ydotdot_rel = delta_ydotdot_new - delta_ydotdot

    delta_z_rel       = delta_z_new - delta_z
    delta_zdot_rel    = delta_zdot_new - delta_zdot
    delta_zdotdot_rel = delta_zdotdot_new - delta_zdotdot

    x       = x0_d + Vx*td + delta_x_rel
    xdot    = Vx + delta_xdot_rel
    xdotdot = delta_xdotdot_rel

    y       = delta_y_rel
    ydot    = delta_ydot_rel
    ydotdot = delta_ydotdot_rel

    z       = delta_z_rel
    zdot    = delta_zdot_rel
    zdotdot = delta_zdotdot_rel

    fx = m * thrust_x(x,zdot,xdotdot,omega,t)
    fy = m * thrust_y(y,ydotdot,omega,t)
    fz = m * thrust_z(z,zdotdot,xdot,omega,t)
    ftot = abs(fx + fy + fz)
    mp = act.RCSpropellant(ftot,dt,Isp)
    m -= mp

    f_array = np.append(f_array,[[fx,fy,fz]],axis=0)
    X_array = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp)
    t_array = np.append(t_array,t)

    delta_x       = delta_x_new
    delta_xdot    = delta_xdot_new
    delta_xdotdot = delta_xdotdot_new

    delta_y       = delta_y_new
    delta_ydot    = delta_ydot_new
    delta_ydotdot = delta_ydotdot_new

    delta_z       = delta_z_new
    delta_zdot    = delta_zdot_new
    delta_zdotdot = delta_zdotdot_new

    t  += dt

#=====================================================================================================================================================================================================================
# Thruster errors
#=====================================================================================================================================================================================================================
def thrust_error(f,cg,angle):
    lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top = thruster_arms(cg)

    T_error_x = np.sin(angle*np.pi/180) * lz_bottom * f
    T_error_y = T_error_x
    T_error_z = np.sin(angle*np.pi/180) * lx_bottom * f

    return T_error_x, T_error_y, T_error_z

def engine_failure():


#Misalignment
cg_orbit = act.z_cg_orbit
angle = 2*np.pi/180

T_error_z1, T_error_y1, T_error_x1 = act.thrust_error(f_array[:,0],cg_orbit,angle)
T_error_z2, T_error_y2, T_error_x2 = act.thrust_error(f_array[:,1],cg_orbit,angle)
T_error_z3, T_error_y3, T_error_x3 = act.thrust_error(f_array[:,2],cg_orbit,angle)

T_error_x_tot = T_error_x1 + T_error_x2 + T_error_x3
T_error_y_tot = T_error_y1 + T_error_y2 + T_error_y3
T_error_z_tot = T_error_z1 + T_error_z2 + T_error_z3
print(T_error_x_tot)
RCS_thrust_error_x = act.RCS_torque_to_thrust(T_error_z_tot,'z',cg_orbit,'error_bottom')
RCS_thrust_error_y = act.RCS_torque_to_thrust(T_error_y_tot,'y',cg_orbit,'error_bottom')
RCS_thrust_error_z = act.RCS_torque_to_thrust(T_error_x_tot,'x',cg_orbit,'error_bottom')
print(RCS_thrust_error_x)
RCS_impulse_error = np.sum(RCS_thrust_error_x*dt) + np.sum(RCS_thrust_error_y*dt) + np.sum(RCS_thrust_error_z*dt)
mp_error           = RCS_impulse_error / (Isp * g_earth)

f_array[:,0] += RCS_thrust_error_x
f_array[:,1] += RCS_thrust_error_y
f_array[:,2] += RCS_thrust_error_z

#Engine failure


#=====================================================================================================================================================================================================================
# Rotations
#=====================================================================================================================================================================================================================
def slew(slew_angle,slew_duration,I):
    slew_acc_duration = 0.05 * slew_duration
    slew_dec_duration = slew_acc_duration
    spin_rate =  slew_angle / slew_duration
    spin_acc  =  spin_rate  / slew_acc_duration
    spin_dec  = -spin_acc
    print('acc',spin_acc)
    RCS_torque = (I * spin_acc)
    return RCS_torque

angle_transfer      = 180 * np.pi / 180
angle_rendezvous    = 180 * np.pi / 180

time_transfer       = 20.
time_rendezvous     = 30.

RCS_torque_transfer   = slew(angle_transfer,time_transfer,Iy)
RCS_thrust_transfer = act.RCS_torque_to_thrust(RCS_torque_transfer,'y',cg_orbit,'normal')
mp_transfer           = 2 * 2 * act.RCSpropellant(RCS_thrust_transfer,time_transfer,Isp)

RCS_torque_rendezvous = slew(angle_rendezvous,time_rendezvous,Ix)
RCS_thrust_rendezvous = act.RCS_torque_to_thrust(RCS_torque_rendezvous,'y',cg_orbit,'normal')
mp_rendezvous         = 2 * act.RCSpropellant(RCS_thrust_rendezvous,time_rendezvous,Isp)

print('rotation',mp_transfer,mp_rendezvous)
print('rotation',RCS_thrust_transfer,mp_rendezvous)

#=====================================================================================================================================================================================================================
# Total thrust and propellant mass
#=====================================================================================================================================================================================================================

#Delta V maneuvers in between phases
thrust_deltaV = thrust_deltaV1 + thrust_deltaV2 + thrust_deltaV3 + thrust_deltaV4 + thrust_deltaV5
mp_deltaV = mp_deltaV1 + mp_deltaV2 + mp_deltaV3 + mp_deltaV4 + mp_deltaV5

#Total
mp_tot     = np.sum(mp) + mp_error

#RCS thrust per engine
RCS_thrust_x    = act.RCS_displacement_to_thrust(f_array[:,0],'z')
RCS_thrust_y    = act.RCS_displacement_to_thrust(f_array[:,1],'y')
RCS_thrust_z    = act.RCS_displacement_to_thrust(f_array[:,2],'x')

print(RCS_thrust_x)

print('Total propellant used: ', mp_tot)
print('Thrust for initial and final delta Vs: ', thrust_deltaV)
print('Max thrust per engine during phases (x,y,z): ', max(RCS_thrust_x), max(RCS_thrust_y), max(RCS_thrust_z))
print('Max error thrust per engine during phases (x,y,z): ', max(RCS_thrust_error_x), max(RCS_thrust_error_y), max(RCS_thrust_error_z))

print(t_array.shape,X_array.shape)
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

    #Propellant mass
    axs[1][0].plot(t_array,mp_array,color="navy")
    axs[1][0].grid(color="gainsboro")
    axs[1][0].set_xlabel("Time [s]")
    axs[1][0].set_ylabel("Propellant mass [kg]")
    plt.show()

    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    fig.suptitle('RCS thrust',fontsize=16)
    #Fx vs time
    axs[0][0].plot(t_array,f_array[:,0],color="navy")
    axs[0][0].grid(color="gainsboro")
    axs[0][0].set_xlabel("Time [s]")
    axs[0][0].set_ylabel("Thrust in X [N]")

    #Fx vs time
    axs[0][1].plot(t_array,f_array[:,1],color="navy")
    axs[0][1].grid(color="gainsboro")
    axs[0][1].set_xlabel("Time [s]")
    axs[0][1].set_ylabel("Thrust in Y [N]")

    #Fx vs time
    axs[0][2].plot(t_array,f_array[:,2],color="navy")
    axs[0][2].grid(color="gainsboro")
    axs[0][2].set_xlabel("Time [s]")
    axs[0][2].set_ylabel("Thrust in Z [N]")
    plt.show()