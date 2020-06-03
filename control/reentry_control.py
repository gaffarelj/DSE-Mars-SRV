import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
from reentry_footprint import flight, time, dt, mean_radius, mars, rho_0
import disturbances as dist
import actuator_properties as act

import numpy as np
from matplotlib import pyplot as plt

#Constants
Iy = act.Iy
cg = act.cg_empty
length = act.length
width = act. width
Isp = act.Isp


#Initial slew values
#Assume spin acceleration/deceleration of 5%, coast time of 90%
slew_duration = 100 #s
slew_coast_duration = 0.90 * slew_duration
slew_acc_duration = 0.05 * slew_duration
alpha = 45 * np.pi / 180 #radians

#flight profile
height = flight[:,3]-mean_radius
velocity = flight[:,0]
gamma = flight[:,1]
rho = mars.density
pitch = gamma + alpha
cd = 2.75
rho = []
for i in height:
    rho.append(mars.density(i))
rho = np.asarray(rho)


#drag disturbance
S, cp = dist.Drag_surface_cp(alpha)
drag = dist.Drag_force(rho,velocity,cd,S)
Td = dist.aerodynamic_disturbance(cp,cg,drag,alpha)


#Torque required for slew during landing
# def slew_landing(alpha,S,cp,rho,velocity,cd,Td,cg,I):
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
#     v0             = velocity[slew_i]
#     Td0            = Td[slew_i]
#
#     spin_rate =  slew_angle_tot / slew_duration
#
# print('thrust per RCS engine:',RCS_thrust_max)
# print('Impulse per RCS engine:', RCS_impulse)
# print('Total propellant needed:', Mp)



def slew_landing(thrust,alpha,S,cp,rho,velocity,cd,Td,cg,I):
    RCS_torque = act.RCS_thrust_to_torque(thrust,"y",length,width,cg)
    slew_angle_tot = np.pi / 180 - alpha
    slew_i         = -int(slew_duration / dt)

    v0             = velocity[slew_i]
    Td0            = Td[slew_i]

    spin_rate =  slew_angle_tot / slew_duration
    spin_acc  =  spin_rate      / slew_acc_duration
    spin_dec  = -spin_acc

    slew_time  = 0
    slew_angle = 0
    impulse    = 0
    while slew_angle < slew_angle_tot:
        S_new, cp_new = dist.Drag_surface_cp(alpha+slew_angle)
        drag = dist.Drag_force(rho[i],velocity[i],cd,S_new)
        Td_new = dist.aerodynamic_disturbance(cp_new,cg,drag_new,alpha+slew_angle)
        dTd = Td_new - Td0
        net_torque = RCS_torque - dTd

        spin_acc    = (net_torque) / I
        spin_rate   = spin_acc * dt
        slew_angle += spin_rate * dt
        slew_time  += dt
        i += 1
        impulse = thrust * dt
    return slew_time, impulse

for i in time:
    print(i)
