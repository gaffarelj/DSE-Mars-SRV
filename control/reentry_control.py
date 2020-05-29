import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
from reentry_footprint import flight, time, dt, mean_radius, mars, rho_0
import disturbances as dist

import numpy as np
from matplotlib import pyplot as plt

#constants
CD = 0.5
cg_full  = 7.7085 #m
cg_empty = 10.0344 #m
alpha = 40 * np.pi / 180 #radians
landing_angle = 90 * np.pi / 180 #radians
slew_duration = 100 #s
#Assume spin acceleration/deceleration of 5%, coast time of 90%
slew_acc_duration = 0.05 * slew_duration
Ix = 2875350.278 #kg/m^2
Iy = 306700.3372 #kg/m^2
Iz = 2875350.278 #kg/m^2

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
print(drag)
#slew
slew_angle_required = landing_angle - pitch
slew_i = -int(slew_duration / dt)

slew_velocity = velocity[slew_i]
slew_angle0 = slew_angle_required[slew_i]
spin_rate = slew_angle0 / slew_duration
spin_acc  = spin_rate / (0.05*slew_duration)
spin_dec  = -spin_acc

Td0 = Td[slew_i]
slew_angle = 0

#ADJUST VELOCITY VECTOR!!!
T = []
for i in range(slew_i,0):
    # print(time[i])
    slew_angle += spin_rate * dt
    # print(slew_angle)
    S_new, cp_new = dist.Drag_surface_cp(alpha+slew_angle)
    drag_new = dist.Drag_force(rho[i],velocity[i],CD,S_new)
    dTd = dist.aerodynamic_disturbance(cp_new,cg_empty,drag_new,alpha+slew_angle)- Td0

    while i < (slew_i + slew_acc_duration*dt):
        T_net = Iy * spin_acc - dTd
    while i > (slew_i + slew_acc_duration*dt) and i < (0 - slew_acc_duration*dt):
        T_net = -dTd
    while i > (0 - slew_acc_duration*dt):
        T_net = Iy * spin_dec - dTd
    T.append(T_net)
T_total = sum(T)
