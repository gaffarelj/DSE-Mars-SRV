import numpy as np
from matplotlib import pyplot as plt

#vehicle constants
length_body =  14.01       #m
length      =  14.01 + 3.6 #m
width = 7.61         #m
Ix = 2875350.278 #kg/m^2
Iy = 306700.3372 #kg/m^2
Iz = 2875350.278 #kg/m^2
cg_full  = 7.7085 #m
cg_empty = 10.0344 #m

#RCS propellant properties
Isp = 150 #mono liquid, N2H2
#Thruster arms
def thruster_arms(cg,length,width):
    lx = width / 2
    ly = lz = length - cg
    return lx, ly, lz

def RCS_torque_to_thrust(T,axis,length,width,cg):
    lx, ly, lz = thruster_arms(cg,length,width)
    # n: number of thrusters to provide torque
    if axis   == "x":
        #bottom RCS
        n_bottom = 2
        #capsule RCS
        n_top = 2
        n = n_bottom + n_top
        f = T / (lx * n_bottom)
    elif axis == "y":
        #bottom RCS
        n_bottom = 3
        #capsule RCS
        n_top = 3
        f = T / (n_top * ly + n_bottom * cg)
        n = n_bottom + n_top
    elif axis == "z":
        #bottom RCS
        n_bottom = 3
        #capsule RCS
        n_top = 3
        f = T / (n_top * lz + n_bottom * cg)
        n = n_bottom + n_top
    return f, n

def RCS_thrust_to_torque(f,axis,length,width,cg):
    lx, ly, lz = thruster_arms(cg,length,width)
    # n: number of thrusters to provide torque
    if axis   == "x":
        #bottom RCS
        n_bottom = 2
        #capsule RCS
        n_top = 2
        n = n_bottom + n_top
        T = f * (lx * n_bottom)
    elif axis == "y":
        #bottom RCS
        n_bottom = 3
        #capsule RCS
        n_top = 3
        T = f * (n_top * ly + n_bottom * cg)
        n = n_bottom + n_top
    elif axis == "z":
        #bottom RCS
        n_bottom = 3
        #capsule RCS
        n_top = 3
        T = f * (n_top * lz + n_bottom * cg)
        n = n_bottom + n_top
    return T

def RCS_displacement_to_thrust(F,axis):
    if axis == "x":
        n_bottom = 4
        n_top    = 4
        n        = n_bottom + n_top
    if axis == "y" or axis == "z":
        n_bottom = 6
        n_top    = 6
        n        = n_bottom + n_top
    f = F / n
    return f, n

def RCSpropellant(Impulse,Isp,n):
    g = 9.81
    Mp = Impulse / (Isp * g) * n
    return Mp
