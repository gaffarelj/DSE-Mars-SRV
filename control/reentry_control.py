import sys
sys.path.append('../astrodynamics')
# import mars_standard_atmosphere as MSA
# from reentry_footprint import flight, time, dt, mean_radius, mars
import disturbances as dist
import actuator_properties as act

import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib import pyplot as plt

margin = 2.
#=====================================================================================================================================================================================================================
#Vehicle onstants
#=====================================================================================================================================================================================================================
Ix = 4105933.742
Iy = 424855.3372
Iz = 4105933.742
cg = act.z_cg_empty
# length = act.length
# width = act. width
Isp = act.Isp_mono
g   = act.g
thrust_levels = np.arange(450,451.85155228408263,0.1)

#=====================================================================================================================================================================================================================
#Flight profile
#=====================================================================================================================================================================================================================
# height = flight[:,3]-mean_radius
# velocity = flight[:,0]

cd = 2.5
t_end = 738.5 #s
q  = 800.  #900. #Pa
alpha = 45 *np.pi / 180
velocity = 500. #m/s
R_mars   = 3389.5*10**3

#=====================================================================================================================================================================================================================
#Disturbance
#=====================================================================================================================================================================================================================
angle = 2. * np.pi /180
omega = velocity / R_mars
S, cp = dist.Drag_surface_cp(angle)
drag = dist.Drag_force(q,cd,S)
Td = dist.aerodynamic_disturbance(cp,cg,drag,angle)
Tgx = dist.gravitygradient_disturbance(Iy,Iz,omega,angle)
Tgy = dist.gravitygradient_disturbance(Ix,Iz,omega,angle)
Tsp= dist.solarpressure_disturbance(angle,cg)
Tm = dist.magnetic_disturbance(R_mars)
T_dist = Tgx + Tsp + Tm


#=====================================================================================================================================================================================================================
#Function to calculate total impulse required for slew maneuver
#=====================================================================================================================================================================================================================
def slew_landing(thrust,alpha0,S,cp,cd,q,T_dist,cg,I,t0,t_end,Isp,g,margin):
    slew_angle_tot = 180 * np.pi / 180 - alpha0
    thrust = thrust
    slew_duration = t_end - t0
    spin_rate_avg = slew_angle_tot / slew_duration
    t             = 0
    dt            = 0.1
    slew_angle    = 0
    spin_rate     = 0
    impulse       = 0
    mp            = 0

    while slew_angle < slew_angle_tot/2 and t < (t_end-t0)/2:
        alpha         = alpha0 + slew_angle
        # if sum(spin_rates)/len(spin_rates) > spin_rate_avg:
        #     thrust = RCS_thrust * 0
        # else:
        #     thrust = RCS_thrust
        RCS_torque    = act.RCS_thrust_to_torque(thrust,"z",cg)

        net_torque    = RCS_torque - T_dist
        # print(T_dist,RCS_torque)
        spin_acc      = (net_torque) / I
        spin_rate    += spin_acc * dt

        slew_angle   += spin_rate * dt
        t            += dt
        impulse      += thrust * dt
        mp           += (thrust * dt) / (Isp*g)

    #Rotation successfully completed or not
    if slew_angle >= slew_angle_tot/2:
        success = 'yes'
    else:
        success = 'no'
    mp = mp * 2
    return impulse,t, mp, success,RCS_torque

#=====================================================================================================================================================================================================================
#All possible rotations with corresponding impulse and time
#=====================================================================================================================================================================================================================

rotation_values = []
for t0 in range(int(t_end-100),int(t_end-10)):

    for thrust in thrust_levels:

            impulse, slew_time, mp, success, RCS_torque = slew_landing(thrust,alpha,S,cp,cd,q,T_dist,cg,Iy,t0,t_end,Isp,g,margin)
            thrust = thrust
            mp     = mp * margin
            if success == 'yes':
                # print([thrust,impulse,slew_time,mp])
                rotation_values.append([thrust,impulse,slew_time,mp,RCS_torque])


rotation_values = np.array(rotation_values)
mp = rotation_values[-1,3]
thrust = rotation_values[-1,0]
t  = rotation_values[-1,2]
torque = rotation_values[-1,-1]

#=====================================================================================================================================================================================================================
#Errors
#=====================================================================================================================================================================================================================
#========================================================================
# Thruster misalignment
#========================================================================
error_angle = 2 * np.pi / 180

T_error_y, T_error_x, T_error_z = act.thrust_error(thrust,cg,error_angle)

RCS_error_x  = act.RCS_torque_to_thrust(T_error_x,'y',cg,'error_bottom')
RCS_error_y  = act.RCS_torque_to_thrust(T_error_y,'x',cg,'error_bottom')
RCS_error_z  = act.RCS_torque_to_thrust(T_error_z,'z',cg,'error_bottom')

RCS_error    = max([RCS_error_x,RCS_error_y,RCS_error_z])
mp_error     = 18 * act.RCSpropellant(RCS_error,t,Isp)

print('==========================================================')
print('REENTRY')
print('THRUST')
print('Total torque (x,y,z)     : ', 0., 0., torque)
print('Thrust per engine (x,y,z): ', 0., 0., thrust)
print('propellant needed: '    , mp)
print('REDUNDANCY')
print('Misalignment torque: ', T_error_x-Tgx-Tsp-Tm,T_error_y-Tgy,T_error_z)
print('Disturbance torque: ', Tgx+Tsp+Tm,Tgy,0)
print('redundancy thrust per engine: ', RCS_error_x,RCS_error_y,RCS_error_z)
print('redundacy propellant per engine: ', mp_error)

#=====================================================================================================================================================================================================================
#Pitch control
#=====================================================================================================================================================================================================================
pitch_moment = np.ones(int(t_end))*14500
pitch_thrust = 2*act.RCS_torque_to_thrust(1000,'z',cg,'normal')
print(pitch_thrust)
pitch_mp     = 0.

# for thrust in pitch_thrust:
#     pitch_mp += 4 * act.RCSpropellant(thrust,1.,Isp)
