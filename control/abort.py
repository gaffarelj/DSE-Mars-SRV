import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
import actuator_properties as act


margin = 2.

p = 101325
A = 0.001
F = p*A

t_burn                  = 100
thrust_level            = 10
Isp_mono                = 140
g                       = 9.80665

capsule_radius_bottom   = 2.3
capsule_radius_top      = 1.4
capsule_height          = 1.2
cg_z                    = 0.3
cg_xy                   = 0.

l_z                     = capsule_height- cg_z
l_z_thrust              = cg_z
l_x = l_y               = capsule_radius_bottom
l_x_thrust = l_y_thrust = capsule_radius_bottom

T_z                     = l_z * F
T_x                     = l_x * F
T_y                     = l_y * F

thrust_z                = T_z / cg_z
thrust_y                = T_y / l_y
thrust_x                = T_x / l_x


def abort_thrust(l,F,Isp,g,l_thrust,thrust_level):
    t  = 0.
    dt = 0.1
    T  = l * F
    mp = 0
    T_net = T
    T_thrust = 2 * thrust_level * l_thrust
    while T_net > 0:
        T_thrust
        T_net = T - T_thrust
        T = T_net
        t += dt
        impulse = 2 * thrust_level* dt * margin
        mp += impulse/(Isp*g)
        print(T_thrust, T_net)

    return t, thrust_level, mp




abort_z = abort_thrust(l_z,F,Isp_mono,g,l_z_thrust,thrust_level)
abort_x = abort_thrust(l_x,F,Isp_mono,g,l_x_thrust,thrust_level)
abort_y = abort_thrust(l_y,F,Isp_mono,g,l_y_thrust,thrust_level)

print(abort_z,abort_x,abort_y)

m_p = abort_z[0] + abort_x[0] + abort_y[0]
