import numpy as np
from matplotlib import pyplot as plt

length = 3.6 + 14.01 #m
width = 7.61         #m
cg_full  = 7.7085 #m
cg_empty = 10.0344 #m

#Thruster arms
def thruster_arms(cg):
    lx = width / 2
    ly = lz = length - cg
    return lx, ly, lz

def RCS_torque_to_thrust(T,axis,length,width,cg):
    lx, ly, lz = thruster_arms(cg)
    # n: number of thrusters to provide torque
    if axis   == "x":
        #bottom RCS
        n_bottom = 2
        f = T / (lx * n_bottom)
    elif axis == "y":
        #bottom RCS
        n_bottom = 3
        #capsule RCS
        n_top = 2

        f = T / (n_top * ly + n_bottom * cg)
    elif axis == "z":
        #bottom RCS
        n_bottom = 2
        #capsule RCS
        n_top = 1
        f = T / (n_top * lz + n_bottom * cg)
    return f
