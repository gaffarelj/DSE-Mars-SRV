import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
from reentry_footprint import flight, time, dt, mean_radius, mars, rho_0
import disturbances as dist
import actuator_properties as act

import numpy as np
from matplotlib import pyplot as plt

#constants
length = 3.6 + 14.01 #m
width = 7.61         #m
Ix = 2875350.278 #kg/m^2
Iy = 306700.3372 #kg/m^2
Iz = 2875350.278 #kg/m^2
CD = 2.75
cg_full  = 7.7085 #m
cg_empty = 10.0344 #m
alpha = 45 * np.pi / 180 #radians

#Initial slew values
#Assume spin acceleration/deceleration of 5%, coast time of 90%
slew_duration = 100 #s
slew_acc_duration = 0.05 * slew_duration


#flight profile
height = flight[:,3]-mean_radius
velocity = flight[:,0]
gamma = flight[:,1]
rho = mars.density
pitch = gamma + alpha
rho = []
for i in height:
    rho.append(mars.density(i))
rho = np.asarray(rho)


#drag disturbance
S, cp = dist.Drag_surface_cp(alpha)
drag = dist.Drag_force(rho,velocity,CD,S)
Td = dist.aerodynamic_disturbance(cp,cg_empty,drag,alpha)


#Torque required for slew during landing
def slew_landing(alpha,S,cp,rho,velocity,CD,Td,cg):
    slew_angle_tot = np.pi / 180 - alpha
    slew_i         = -int(slew_duration / dt)

    v0             = velocity[slew_i]
    Td0            = Td[slew_i]

    spin_rate =  slew_angle_tot / slew_duration
    spin_acc  =  spin_rate      / slew_acc_duration
    spin_dec  = -spin_acc


    Torque = []
    slew_angle = 0
    for i in range(slew_i,0):
        slew_angle += spin_rate * dt
        # print(slew_angle)
        S_new, cp_new = dist.Drag_surface_cp(alpha+slew_angle)
        drag_new      = dist.Drag_force(rho[i],velocity[i],CD,S_new)
        Td_new        = dist.aerodynamic_disturbance(cp_new,cg,drag_new,alpha+slew_angle)
        dTd           = Td_new - Td0
        # print(Td_new)
        if i < (slew_i + slew_acc_duration/dt):
            Torque_net = Iy * spin_acc - dTd
        elif i > (slew_i + slew_acc_duration/dt) and i < (0 - slew_acc_duration/dt):
            Torque_net = -dTd
        elif i > (0 - slew_acc_duration/dt):
            Torque_net = Iy * spin_dec - dTd
        Torque.append(Torque_net)
    Torque = np.array(Torque)
    Impulse = np.sum(Torque*dt)
    return Torque


slew_torque_landing = slew_landing(alpha,S,cp,rho,velocity,CD,Td,cg_empty)
print(slew_torque_landing)

RCS_thrust  = act.RCS_torque_to_thrust(slew_torque_landing,"y",length,width,cg_empty)
RCS_impulse = np.sum(RCS_thrust * dt)
print(RCS_thrust,RCS_impulse)



def slew_landing(disturbance, alpha,S,cp,rho,velocity,CD,Td,cg):
    slew_angle_tot = np.pi / 180 - alpha
    slew_i         = -int(slew_duration / dt)

    spin_rate =  slew_angle_tot / slew_duration
    spin_acc  =  spin_rate      / slew_acc_duration
    spin_dec  = -spin_acc


    Torque = []
    slew_angle = 0
    for i in range(slew_i,0):
        slew_angle += spin_rate * dt
        # print(slew_angle)
        if disturbance != False:
            S_new, cp_new = dist.Drag_surface_cp(alpha+slew_angle)
            drag_new      = dist.Drag_force(rho[i],velocity[i],CD,S_new)
            Td_new        = dist.aerodynamic_disturbance(cp_new,cg,drag_new,alpha+slew_angle)
            dTd           = Td_new - Td0
        # print(Td_new)
        if i < (slew_i + slew_acc_duration/dt):
            Torque_net = Iy * spin_acc - dTd
        elif i > (slew_i + slew_acc_duration/dt) and i < (0 - slew_acc_duration/dt):
            Torque_net = -dTd
        elif i > (0 - slew_acc_duration/dt):
            Torque_net = Iy * spin_dec - dTd
        Torque.append(Torque_net)
    Torque = np.array(Torque)
    Impulse = np.sum(Torque*dt)
    return Torque


slew_torque_landing = slew_landing(alpha,S,cp,rho,velocity,CD,Td,cg_empty)
print(slew_torque_landing)

RCS_thrust  = act.RCS_torque_to_thrust(slew_torque_landing,"y",length,width,cg_empty)
RCS_impulse = np.sum(RCS_thrust * dt)
print(RCS_thrust,RCS_impulse)
