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
Isp = 140 #bi liquid, LCH4
Isp_mono = 140 #mono liquid H2O2
g   = 9.80665

#engines
nx = 6
ny = 6
nz = 4

def thruster_arms(z_cg):
    #vehicle constants
    length_body             =  14.01
    length                  =  14.01 + 3.6
    body_radius             = 7.61 / 2
    capsule_radius_bottom   = 2.3
    capsule_radius_top      = 1.4

    x_cg                    = 0.
    y_cg                    = 0.
    y_cg_full               = 7.7085 #m
    y_cg_empty              = 10.0344 #m
    y_cg_orbit              = 6.62390

    x_body_side             = body_radius
    x_capsule_bottomside    = capsule_radius_bottom
    x_capsule_topside       = capsule_radius_top

    z_body_side             = body_radius
    z_capsule_bottomside    = capsule_radius_bottom
    z_capsule_topside       = capsule_radius_top

    y_body_top              = length_body
    y_capsule_bottom        = y_body_top   + 1.2
    y_capsule_top           = y_capsule_bottom + 2.7

    lx_bottom = x_body_side
    lx_top    = x_capsule_bottomside

    lz_bottom = z_body_side
    lz_top    = z_capsule_bottomside

    ly_bottom = y_cg - 1.
    ly_top    = y_cg_orbit * 2 - y_cg

    return lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top

def RCS_torque_to_thrust(T,axis,cg,scenario):
    lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top = thruster_arms(cg)
    # n: number of thrusters to provide torque

    if axis   == 'x' or axis == 'z':
        if scenario == 'normal':
            n_bottom = 2
            n_top = 2
            thrust = T / (ly_bottom * n_bottom + ly_top * n_top)
        elif scenario == 'error_bottom' or scenario == 'error_top':
            n = 1
            thrust = T / (lz_bottom * n)
        elif scenario == 'failure':
            n_bottom = 1
            n_top = 1
            thrust = T / (lz_bottom * n_bottom + lz_top * n_top)

    elif axis == "y":
        if scenario == 'normal':
            n_bottom = 2
            n_top = 2
            thrust = T / (n_top * ly_top + n_bottom * ly_bottom)
        elif scenario == 'error_bottom':
            n = 1
            thrust = T / (lx_bottom)
        elif scenario == 'error_top':
            n = 1
            thrust = T / (lx_top)
        elif scenario == 'failure':
            n_bottom = 1
            n_top = 1
            thrust = T / (lz_bottom * n_bottom + lz_top * n_top)

    return thrust


def RCS_displacement_to_thrust(F,axis,scenario):
    if axis == "x" or axis == 'z':
        if scenario == 'normal':
            n_bottom = 1
            n_top    = 1
            n        = n_bottom + n_top
        elif scenario == 'failure':
            n_bottom = 1
            n_top    = 1
            n        = n_bottom + n_top

    if axis == "y":
        if scenario == 'normal':
            n_bottom = 1000
            n_top    = 0
            n        = n_bottom + n_top
        if scenario == 'failure':
            n_bottom = 2
            n_top    = 0
            n        = n_bottom + n_top
    f = F / n
    return f

def RCS_thrust_to_torque(f,axis,cg):
    lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top = thruster_arms(cg)
    if axis   == 'x' or axis == 'z':
        n_bottom = 2
        n_top = 2
        T = f * (lz_bottom * n_bottom + lz_top * n_top)

    elif axis == 'y':
        #bottom RCS
        n_bottom = 2
        n_top = 2
        T = f * (n_top * ly_top + n_bottom * ly_bottom)
    return T

def slew(thrust,tburn,slew_angle,I):
    torque = RCS_thrust_to_torque(thrust,'z','normal')
    spin_acc = torque / I
    spin_rate = spin_acc * tburn
    slew_time = slew_angle / spin_rate

    return slew_time,torque

def thrust_error(f,cg,angle):
    lx_bottom, lx_top, ly_bottom, ly_top, lz_bottom, lz_top = thruster_arms(cg)

    T_error_x = np.sin(angle*np.pi/180) * ly_bottom * f
    T_error_z = T_error_x
    T_error_y = np.sin(angle*np.pi/180) * lx_bottom * f

    return T_error_x, T_error_y, T_error_z


def RCSpropellant(f,t,Isp):
    g = 9.80665
    impulse = f * t
    Mp = impulse / (Isp * g)
    return Mp
