import numpy as np
from matplotlib import pyplot as plt

#vehicle constants
length_body             =  14.01
length                  =  14.01 + 3.6
body_radius             = 7.61 / 2
capsule_radius_bottom   = 2.3
capsule_radius_top      = 1.4

x_cg                    = 0.
z_cg                    = 0.
z_cg_full               = 7.7085 #m
z_cg_empty              = 10.0344 #m
z_cg_orbit              = 6.62390

x_body_side             = body_radius
x_capsule_bottomside    = capsule_radius_bottom
x_capsule_topside       = capsule_radius_top

y_bdoy_side             = body_radius
y_capsule_bottomside    = capsule_radius_bottom
y_capsule_topside       = capsule_radius_top

z_body_top              = length_body
z_capsule_bottom        = z_body_top   + 1.2
z_capsule_top           = z_capsule_bottom + 2.7


Ix = 2875350.278 #kg/m^2
Iz = 306700.3372 #kg/m^2
Iy = 2875350.278 #kg/m^2

#RCS propellant properties
Isp = 390 #mono liquid, N2H2

def thruster_arms(z_cg):
    #vehicle constants
    length_body             =  14.01
    length                  =  14.01 + 3.6
    body_radius             = 7.61 / 2
    capsule_radius_bottom   = 2.3
    capsule_radius_top      = 1.4

    x_cg                    = 0.
    y_cg                    = 0.
    z_cg_full               = 7.7085 #m
    z_cg_empty              = 10.0344 #m
    z_cg_orbit              = 6.62390

    x_body_side             = body_radius
    x_capsule_bottomside    = capsule_radius_bottom
    x_capsule_topside       = capsule_radius_top

    y_body_side             = body_radius
    y_capsule_bottomside    = capsule_radius_bottom
    y_capsule_topside       = capsule_radius_top

    z_body_top              = length_body
    z_capsule_bottom        = z_body_top   + 1.2
    z_capsule_top           = z_capsule_bottom + 2.7

    lx_bottom = x_body_side
    lx_top    = x_capsule_bottomside

    ly_bottom = y_body_side
    ly_top    = y_capsule_bottomside

    lz_bottom = z_cg
    lz_top    = z_cg_orbit * 2 - z_cg

    return lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top

def RCS_torque_to_thrust(T,axis,cg,scenario):
    lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top = thruster_arms(cg)
    # n: number of thrusters to provide torque

    if axis   == 'x' or axis == 'y':
        if scenario == 'normal':
            n_bottom = 3
            n_top = 3
            thrust = T / (lz_bottom * n_bottom + lz_top * n_top)
        if scenario == 'error_bottom' or scenario == 'error_top':
            n = 1
            thrust = T / (lz_bottom * n)

    elif axis == "z":
        if scenario == 'normal':
            n_bottom = 2
            n_top = 2
            thrust = T / (n_top * lx_top + n_bottom * lx_bottom)
        elif scenario == 'error_bottom':
            n = 1
            thrust = T / (lx_bottom)
        elif scenario == 'error_top':
            n = 1
            thrust = T / (lx_top)

    return thrust


def RCS_displacement_to_thrust(F,axis):
    if axis == "x" or axis == 'y':
        n_bottom = 3
        n_top    = 3
        n        = n_bottom + n_top
    if axis == "y" or axis == "z":
        n_bottom = 4
        n_top    = 0
        n        = n_bottom + n_top
    f = F / n
    return f

def RCS_thrust_to_torque(f,axis,cg):
    lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top = thruster_arms(cg)
    if axis   == 'x' or axis == 'y':
        n_bottom = 3
        n_top = 3
        T = f * (lz_bottom * n_bottom + lz_top * n_top)

    elif axis == "z":
        #bottom RCS
        n_bottom = 2
        n_top = 2
        T = f * (n_top * lx_top + n_bottom * lx_bottom)

    return T



def RCSpropellant(f,t,Isp):
    g = 9.80665
    impulse = f * t
    Mp = impulse / (Isp * g)
    return Mp
