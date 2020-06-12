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

mu     = 0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
R      = 3389.5*10**3                  #[m] volumetric mean radius
h_node = 500*10**3
h_phasing=609.74*10**3
g_earth= 9.81             #[m]

period =  2*np.pi* np.sqrt((R+h_node) ** 3 / mu)
omega  = (2 * np.pi)/period
time_hohmann = np.pi*np.sqrt(((h_node+h_phasing+2*R)/2)**3*1/mu)
V      = np.sqrt((R+h_node) / mu)

#=====================================================================================================================================================================================================================
#Vehicle properties
#=====================================================================================================================================================================================================================

#Implement Mass of Charon at t0!
m0   =53248.16053543461         #previous value was: 100000
Isp  = act.Isp
Isp_mono = act.Isp         #previous value was: 390
OF_ratio = 3.8
Ix = 4069666.12 #kg/m^2
Iy = 276844.4638 #kg/m^2
Iz = 4069666.12 #kg/m^2
cg_orbit = act.z_cg_orbit

#=====================================================================================================================================================================================================================
#Vbar approach properties
#=====================================================================================================================================================================================================================
margin = 2.

#Proximity operations A:
x0_A = -1000         #m
x1_A = -250          #m
t_A  =  period   #s
Vx_A   = (x1_A - x0_A) / t_A #m/s
deltaV_A_0 = deltaV_A_1 = Vx_A * margin

#Proximity operations B:
x0_B = x1_A         #m
x1_B = -30           #m
t_B  = 60 * 60      #s
Vx_B  = (x1_B - x0_B) / t_B #m/s
deltaV_B_0 = deltaV_B_1 = Vx_B * margin
#Docking:
x0_d = x1_B         #m
x1_d = -3            #m
t_d  = (5*((x1_d - x0_d)/10)) * 60       #m
Vx_d  = (x1_d - x0_d) / t_d #m/s
deltaV_d_0 = deltaV_d_1 = Vx_d * margin

y0 = z0 = 0

deltaV_tot = deltaV_A_0 + deltaV_A_1 + deltaV_B_0 + deltaV_B_1 + deltaV_d_0 + deltaV_d_1

print(2*omega*Vx_A*m0,2*omega*Vx_B*m0,2*omega*Vx_d*m0)
#=====================================================================================================================================================================================================================
#Navigation measurement errors
#=====================================================================================================================================================================================================================
error_xm = 10.
error_ym = 10.
error_zm = 10.

error_xdot = 1.
error_ydot = 1.
error_zdot = 1.

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

plotting=True      #Do you wanna plot? no=False

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
tburn = 6.
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
# print('initial xyz: ', x,y,z)
# print('initial delta xyz: ', delta_x,delta_y,delta_z)

fx0 = m * thrust_x(x,zdot,xdotdot,omega,t)
fy0 = m * thrust_y(y,ydotdot,omega,t)
fz0 = m * thrust_z(z,zdotdot,xdot,omega,t)

f_array   = np.array([[fx0,fy0,fz0]])
X_array   = np.array([[x,y,z]])
t_array   = np.array(t)
mp_array  = np.array(mp)
mp_biprop_array = np.array(mp)

#DeltaV maneuver 1
# t += tburn
thrust_deltaV1, mp_deltaV1 = vac_thrust(deltaV_A_0,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV1
f_array = np.append(f_array,[[thrust_deltaV1,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV1)
mp_biprop_array = np.append(mp_biprop_array,mp_deltaV1)
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

    fx = m * thrust_x(x,zdot,xdotdot,omega,t) * margin
    fy = m * thrust_y(y,ydotdot,omega,t) * margin
    fz = m * thrust_z(z,zdotdot,xdot,omega,t) * margin
    ftot = abs(fx + fy + fz)
    mp = act.RCSpropellant(ftot,dt,Isp_mono)
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
# t += tburn
thrust_deltaV2, mp_deltaV2 = vac_thrust(deltaV_A_1,Isp,m,tburn,De=0,pe=0)
m -= mp_deltaV2
f_array = np.append(f_array,[[thrust_deltaV2,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV2)
mp_biprop_array = np.append(mp_biprop_array,mp_deltaV2)
t_array  = np.append(t_array,t)

#Delta V maneuver 3
# t += tburn
thrust_deltaV3, mp_deltaV3 = vac_thrust(deltaV_B_0,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV3
f_array = np.append(f_array,[[thrust_deltaV3,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV3)
mp_biprop_array = np.append(mp_biprop_array,mp_deltaV3)
t_array  = np.append(t_array,t)

#Proximity operations B
while x >= x0_B-1 and x < x1_B:
    t  += dt
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

    fx = m * thrust_x(x,zdot,xdotdot,omega,t) * margin
    fy = m * thrust_y(y,ydotdot,omega,t) * margin
    fz = m * thrust_z(z,zdotdot,xdot,omega,t) * margin
    ftot = abs(fx + fy + fz)
    mp = act.RCSpropellant(ftot,dt,Isp_mono)
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



#Delta V maneuver 4
# t += tburn
thrust_deltaV4, mp_deltaV4 = vac_thrust(deltaV_B_1,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV4
f_array = np.append(f_array,[[thrust_deltaV4,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV4)
mp_biprop_array = np.append(mp_biprop_array,mp_deltaV4)
t_array  = np.append(t_array,t)

#Delta V maneuver 5
# t += tburn
thrust_deltaV5, mp_deltaV5 = vac_thrust(deltaV_d_0,Isp,m0,tburn,De=0,pe=0)
m -= mp_deltaV5
f_array = np.append(f_array,[[thrust_deltaV5,0,0]],axis=0)
X_array  = np.append(X_array,[[x,y,z]],axis=0)
mp_array = np.append(mp_array,mp_deltaV5)
mp_biprop_array = np.append(mp_biprop_array,mp_deltaV5)
t_array  = np.append(t_array,t)

while x >= x0_d-1 and x < x1_d:
    Vx  = Vx_d
    td += dt
    t  += dt

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

    fx = m * thrust_x(x,zdot,xdotdot,omega,t) * margin
    fy = m * thrust_y(y,ydotdot,omega,t) * margin
    fz = m * thrust_z(z,zdotdot,xdot,omega,t) * margin
    ftot = abs(fx + fy + fz)
    mp = act.RCSpropellant(ftot,dt,Isp_mono)
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



#=================================================================================================================================================
# Total thrust and propellant mass
#=================================================================================================================================================

#Delta V maneuvers in between phases
thrust_deltaV = [thrust_deltaV1, thrust_deltaV2, thrust_deltaV3, thrust_deltaV4, thrust_deltaV5]
mp_deltaV = mp_deltaV1 + mp_deltaV2 + mp_deltaV3 + mp_deltaV4 + mp_deltaV5

#Total
mp_tot     = np.sum(mp_array)

#RCS thrust per engine
RCS_thrust_x    = act.RCS_displacement_to_thrust(f_array[:,0],'y','normal')
RCS_thrust_y    = act.RCS_displacement_to_thrust(f_array[:,1],'x','normal')
RCS_thrust_z    = act.RCS_displacement_to_thrust(f_array[:,2],'z','normal')

#=================================================================================================================================================
# Errors
#=================================================================================================================================================
#========================================================================
# Thruster misalignment
#========================================================================
#Misalignment
angle = 2*np.pi/180

T_error_y1, T_error_x1, T_error_z1 = act.thrust_error(max(f_array[:,0]),cg_orbit,angle)
T_error_y2, T_error_x2, T_error_z2 = act.thrust_error(max(f_array[:,1]),cg_orbit,angle)
T_error_y3, T_error_x3, T_error_z3 = act.thrust_error(max(f_array[:,2]),cg_orbit,angle)
T_error_x = T_error_x1 + T_error_x2 + T_error_x3
T_error_y = T_error_y1 + T_error_y2 + T_error_y3
T_error_z = T_error_z1 + T_error_z2 + T_error_z3


#========================================================================
# Disturbances
#========================================================================
theta = 2. * np.pi / 180

Tgx  = dist.gravitygradient_disturbance(Iy,Iz,omega,theta)
Tgy  = dist.gravitygradient_disturbance(Ix,Iz,omega,theta)
Tsp = dist.solarpressure_disturbance(theta,cg_orbit)
Tm  = dist.magnetic_disturbance(R)

#========================================================================
# Total
#========================================================================
T_error_x_tot = T_error_x + Tgx + Tsp + Tm
T_error_y_tot = T_error_y + Tgy
T_error_z_tot = T_error_z

RCS_thrust_error_x = act.RCS_torque_to_thrust(T_error_y_tot,'y',cg_orbit,'error_bottom') * margin
RCS_thrust_error_y = act.RCS_torque_to_thrust(T_error_x_tot,'x',cg_orbit,'error_bottom') * margin
RCS_thrust_error_z = act.RCS_torque_to_thrust(T_error_z_tot,'z',cg_orbit,'error_bottom') * margin

RCS_impulse_error = (np.sum(RCS_thrust_error_x*t) + np.sum(RCS_thrust_error_y*t) + np.sum(RCS_thrust_error_z*t))
mp_error           = RCS_impulse_error / (Isp_mono * g_earth)

#=================================================================================================================================================
# Rotations
#=================================================================================================================================================


angle_transfer      = 180 * np.pi / 180
angle_rendezvous    = 180 * np.pi / 180

time_transfer       = time_hohmann
time_rendezvous     = time_hohmann

RCS_thrust_rotations  = 451.85155228408263
tburn                 = 1.


transfer_time         = act.slew(RCS_thrust_rotations,tburn,angle_transfer,Iy)
mp_rotations          = 4 * 3 * act.RCSpropellant(RCS_thrust_rotations,tburn,Isp_mono)
print('Rotation time limit: ', time_hohmann)
print('Rotation duration  :',transfer_time)
T_error_rot_y, T_error_rot_x, T_error_rot_z = act.thrust_error(RCS_thrust_rotations,cg_orbit,angle)
T_error_rot_x += Tgx + Tsp + Tm
RCS_error_rot_x  = act.RCS_torque_to_thrust(T_error_rot_x,'z',cg_orbit,'error_bottom')
RCS_error_rot_y  = act.RCS_torque_to_thrust(T_error_rot_y,'z',cg_orbit,'error_bottom')
RCS_error_rot_z  = act.RCS_torque_to_thrust(T_error_rot_z,'z',cg_orbit,'error_bottom')
RCS_error_rot    = margin * max(RCS_error_rot_x, RCS_error_rot_y, RCS_error_rot_z)
mp_error_rot     = 18 * act.RCSpropellant(RCS_error_rot,tburn,Isp_mono)


#=================================================================================================================================================
# Final thrust and propellant values
#=================================================================================================================================================
mp_biprop_fuel= np.sum(mp_biprop_array)
mp_biprop_ox  = np.sum(mp_biprop_array) * OF_ratio
mp_total = mp_tot + mp_error + mp_rotations + mp_error_rot
print('==========================================================')
print('==========================================================')
print('APPROACH PHASE')
print('==========================================================')
print('PROPELLANT')
print('Total propellant used: ', mp_total)
print('Redundancy propellant: ', mp_error)
print('==========================================================')
print('Thrust for initial and final delta Vs: ', thrust_deltaV)
print('==========================================================')
print('THRUST')
print('Max thrust total      (x,y,z): ', 4*max(RCS_thrust_x), 2*max(RCS_thrust_y), 2*max(RCS_thrust_z))
print('Max thrust per engine (x,y,z): ', max(RCS_thrust_x), max(RCS_thrust_y), max(RCS_thrust_z))
print('Min thrust total      (x,y,z): ', 4*min(RCS_thrust_x), 2*min(RCS_thrust_y), 2*min(RCS_thrust_z))
print('Min thrust per engine (x,y,z): ', min(RCS_thrust_x), min(RCS_thrust_y), min(RCS_thrust_z))
print('MISALIGNMENT REDUNDANCY')
print('Misalignment torque (x,y,z): '              , T_error_x, T_error_y, T_error_z)
print('Disturbance torque  (x,y,z): '              , Tgx+Tsp+Tm, Tgy, 0.)
print('Redundancy thrust (misalignment in x,y,z): ', RCS_thrust_error_x, RCS_thrust_error_y, RCS_thrust_error_z)
# print('Redundancy thrust (engine failure in x,y,z): ', max(RCS_thrust_failure_x), max(RCS_thrust_failure_y), max(RCS_thrust_failure_z))
print('==========================================================')
print('==========================================================')
print('ROTATIONS')
print('Propellant used (transfer, go-nogo): '      , mp_rotations)
print('Thrust                             : '      , RCS_thrust_rotations)
print('MISALIGNMENT REDUNDANCY')
print('Misalignment torque (x,y,z): '              , T_error_rot_x+2*(Tgx+Tsp+Tm), T_error_rot_y, T_error_rot_z)
print('Disturbance torque  (x,y,z): '              , Tgx+Tsp+Tm, Tgy, 0.)
print('Redundancy thrust (transfer, go-nogo): '    , RCS_error_rot_x + RCS_error_rot_y, RCS_error_rot_z)
print('Redundancy propellant (transfer, go-nogo): ', mp_error_rot)
print('==========================================================')
print('==========================================================')

#=====================================================================================================================================================================================================================
# Plotting
#=====================================================================================================================================================================================================================
if plotting:
    fig, axs = plt.subplots(2, 3, constrained_layout=True)
    fig.suptitle('Relative position wrt target',fontsize=16)
    #Fx vs time
    axs[0][0].plot(t_array[2:],X_array[2:,0],color="navy")
    axs[0][0].grid(color="gainsboro")
    axs[0][0].set_xlabel("Time [s]")
    axs[0][0].set_ylabel("X position [m]")

    #Fx vs time
    axs[0][1].plot(t_array[2:],X_array[2:,1],color="navy")
    axs[0][1].grid(color="gainsboro")
    axs[0][1].set_xlabel("Time [s]")
    axs[0][1].set_ylabel("Y position [m]")

    #Fx vs time
    axs[0][2].plot(t_array[2:],X_array[2:,2],color="navy")
    axs[0][2].grid(color="gainsboro")
    axs[0][2].set_xlabel("Time [s]")
    axs[0][2].set_ylabel("Z postion [m]")
    plt.show()

    fig, axs = plt.subplots(2, 2, constrained_layout=True)
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
    axs[1][0].plot(t_array,f_array[:,2],color="navy")
    axs[1][0].grid(color="gainsboro")
    axs[1][0].set_xlabel("Time [s]")
    axs[1][0].set_ylabel("Thrust in Z [N]")

    #Propellant mass
    axs[1][1].plot(t_array,mp_array,color="navy")
    axs[1][1].grid(color="gainsboro")
    axs[1][1].set_xlabel("Time [s]")
    axs[1][1].set_ylabel("Propellant mass [kg]")
    plt.show()
