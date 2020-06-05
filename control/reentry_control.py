import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
from reentry_footprint import flight, time, dt, mean_radius, mars
import disturbances as dist
import actuator_properties as act

import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib import pyplot as plt

#=====================================================================================================================================================================================================================
#Vehicle onstants
#=====================================================================================================================================================================================================================
Iy = act.Iy
cg = act.cg_empty
# length = act.length
# width = act. width
Isp = act.Isp
g   = act.g
thrust_levels = np.arange(100,1400,200)

#=====================================================================================================================================================================================================================
#Flight profile
#=====================================================================================================================================================================================================================
# height = flight[:,3]-mean_radius
velocity = flight[:,0]

cd = 2.75
t_end = 738.5 #s
q  = 800.  #900. #Pa
alpha = 45 *np.pi / 180
velocity = 500. #m/s

#=====================================================================================================================================================================================================================
#Drag disturbance
#=====================================================================================================================================================================================================================

S, cp = dist.Drag_surface_cp(alpha)
drag = dist.Drag_force(rho,velocity,cd,S)
Td = dist.aerodynamic_disturbance(cp,cg,drag,alpha)

# print(Td)
# # Torque required for slew during landing
# def slew_landing(alpha,S,cp,q,cd,Td,cg,I):
#     slew_angle_tot = np.pi / 180 - alpha
#     slew_i         = -int(slew_duration / dt)
#
#     v0             = velocity[slew_i]
#     Td0            = Td[slew_i]
#
#     spin_rate =  slew_angle_tot / slew_duration
#     spin_acc  =  spin_rate      / slew_acc_duration
#     spin_dec  = -spin_acc
#
#
#     RCS_torque = []
#     slew_angle = 0
#     for i in range(slew_i,0):
#         slew_angle += spin_rate * dt
#         # print(slew_angle)
#         S_new, cp_new = dist.Drag_surface_cp(alpha+slew_angle)
#         drag_new      = dist.Drag_force(rho[i],velocity[i],cd,S_new)
#         Td_new        = dist.aerodynamic_disturbance(cp_new,cg,drag_new,alpha+slew_angle)
#         dTd           = Td_new - Td0
#         # print(Td_new)
#         if i < (slew_i + slew_acc_duration/dt):
#             Torque_net = I * spin_acc - dTd
#         elif i > (slew_i + slew_acc_duration/dt) and i < (0 - slew_acc_duration/dt):
#             Torque_net = -dTd
#         elif i > (0 - slew_acc_duration/dt):
#             Torque_net = I * spin_dec - dTd
#         RCS_torque.append(Torque_net)
#     RCS_torque = np.array(RCS_torque)
#     return RCS_torque
#
# RCS_torque = slew_landing(alpha,S,cp,rho,velocity,cd,Td,cg,Iy)
# RCS_thrust, number_thrusters  = act.RCS_torque_to_thrust(RCS_torque,"y",length,width,cg)
# RCS_impulse = np.sum(RCS_thrust * dt)
# Mp = act.RCSpropellant(RCS_impulse,Isp,number_thrusters)
# RCS_thrust_max = np.max(RCS_thrust)
#
#
# print('thrust per RCS engine:',RCS_thrust_max)
# print('Impulse per RCS engine:', RCS_impulse)
# print('Total propellant needed:', Mp)


#=====================================================================================================================================================================================================================
#Function to calculate total impulse required for slew maneuver
#=====================================================================================================================================================================================================================

def slew_landing(RCS_thrust,alpha,S,cp,cd,q,Td0,cg,I,t0):
    slew_angle_tot = 180 * np.pi / 180 - alpha

    slew_duration = t1 - t0
    spin_rate_avg = slew_angle_tot / slew_duration
    t             = 0
    dt            = 0.01
    slew_angle    = 0
    impulse       = 0
    spin_rates    = []
    spin_rates
    while t < t1 and slew_angle:
        if sum(spin_rates)/len(spin_rates) > spin_rate_avg:
            thrust = RCS_thrust * 0
        else:
            thrust = RCS_thrust

        RCS_torque    = act.RCS_thrust_to_torque(thrust,"y",cg)

        S_new, cp_new = dist.Drag_surface_cp(alpha)
        drag_new      = dist.Drag_force(q,cd,S_new)
        Td_new        = dist.aerodynamic_disturbance(cp_new,cg,drag_new,alpha)
        dTd           = Td_new - Td0

        net_torque    = RCS_torque - dTd
        spin_acc      = (net_torque) / I
        spin_rate     = spin_acc * dt


        slew_angle   += spin_rate * dt
        alpha        += slew_angle
        t            += dt
        impulse      += thrust * dt
        spin_rates.append(spin_rate)

    return impulse,t

#=====================================================================================================================================================================================================================
#All possible rotations with corresponding impulse and time
#=====================================================================================================================================================================================================================

rotation_values = []

for t0 in range(t_end-100,t_end-20):

    for thrust in thrust_levels:

            impulse, slew_time = slew_landing(thrust,alpha,S,cp,cd,q,cg,Iy,t0)
            mp       = impulse / (Isp * g)
            rotation_values.append(impulse,thrust,slew_time,mp)

print(rotation_values)


# tot_times_and_impulses = []
# for i in range(len(time)):
#     t0 = time[i]
#     Td0 = Td[i]
#     V0 = velocity[i]
#     times_and_impulses = []
#
#     for j in thrust_levels:
#         slew_time, impulse = slew_landing(j,alpha,S,cp,cd,rho,V0,Td0,cg,Iy,i)
#         if slew_time > time[-1] - t0:
#             t1 = 0
#         else:
#             t1 = t0 + time[int(slew_time/dt)]
#         times_and_impulses.append((t0,t1,slew_time,impulse,j))
#
#     index_min_impulse = np.argmin(np.array(times_and_impulses)[:,-2])
#     min_impulse = times_and_impulses[index_min_impulse]
#     tot_times_and_impulses.append(min_impulse)
#
# index_tot_min_impulse = np.argmin(np.array(tot_times_and_impulses)[:,-2])
# tot_min_impulse = tot_times_and_impulses[index_tot_min_impulse]
# print(tot_min_impulse)
