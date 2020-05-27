#Made by the aerodynamicist

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import matplotlib.pyplot as plt

#I'm genuinely sorry for how awful this looks. I thought this way would be the neatest, I was wrong.
SS_alpha_spacing = np.linspace(-10,25,8) #measurement every 5 deg
SS_Cl_M_0025 = np.array([-0.50,-0.30,-0.05,0.18,0.40,0.63,0.93,1.10])
SS_Cl_M_0060 = np.array([-0.52,-0.30,-0.05,0.18,0.42,0.67,0.92,1.05])
SS_Cl_M_0080 = np.array([-0.58,-0.35,-0.05,0.19,0.44,0.67,0.90,0.92])
SS_Cl_M_0092 = np.array([-0.60,-0.36,-0.05,0.21,0.48,0.70,0.92,0.96])
SS_Cl_M_0098 = np.array([-0.64,-0.37,-0.05,0.23,0.52,0.78,1.00,1.08])
SS_Cl_M_0110 = np.array([-0.64,-0.38,-0.05,0.24,0.52,0.78,1.00,1.10])
SS_Cl_M_0150 = np.array([-0.48,-0.24,-0.05,0.21,0.42,0.62,0.83,0.99])
SS_Cl_M_0200 = np.array([-0.37,-0.20,-0.05,0.15,0.31,0.50,0.68,0.83])
SS_Cl_M_0300 = np.array([-0.28,-0.18,-0.05,0.08,0.22,0.37,0.50,0.67])
SS_Cl_M_0400 = np.array([-0.22,-0.16,-0.05,0.06,0.19,0.32,0.46,0.62])
SS_Cl_M_0500 = np.array([-0.20,-0.10,-0.05,0.04,0.16,0.28,0.42,0.58])
SS_Cl_M_0800 = np.array([-0.20,-0.10,-0.05,0.01,0.11,0.22,0.37,0.52])
SS_Cl_M_1000 = np.array([-0.20,-0.10,-0.05,0.01,0.10,0.21,0.36,0.50])
SS_Cl_M_1500 = np.array([-0.20,-0.10,-0.05,0.01,0.10,0.21,0.36,0.50])
SS_Cl_M_2000 = np.array([-0.20,-0.10,-0.05,0.01,0.10,0.21,0.36,0.50])

SS_Cd_M_0025 = np.array([0.11,0.08,0.06,0.06,0.10,0.15,0.29,0.45])
SS_Cd_M_0060 = np.array([0.13,0.08,0.06,0.06,0.10,0.19,0.36,0.50])
SS_Cd_M_0080 = np.array([0.16,0.10,0.07,0.07,0.13,0.24,0.39,0.50])
SS_Cd_M_0092 = np.array([0.21,0.14,0.10,0.11,0.17,0.28,0.44,0.56])
SS_Cd_M_0098 = np.array([0.25,0.17,0.14,0.15,0.23,0.35,0.50,0.65])
SS_Cd_M_0110 = np.array([0.28,0.19,0.16,0.18,0.24,0.36,0.52,0.68])
SS_Cd_M_0150 = np.array([0.26,0.19,0.16,0.17,0.22,0.31,0.44,0.59])
SS_Cd_M_0200 = np.array([0.21,0.16,0.14,0.15,0.18,0.25,0.36,0.50])
SS_Cd_M_0300 = np.array([0.18,0.14,0.11,0.11,0.13,0.18,0.27,0.40])
SS_Cd_M_0400 = np.array([0.16,0.12,0.09,0.09,0.11,0.16,0.25,0.36])
SS_Cd_M_0500 = np.array([0.15,0.12,0.09,0.08,0.10,0.15,0.22,0.35])
SS_Cd_M_0800 = np.array([0.14,0.11,0.08,0.07,0.08,0.13,0.19,0.32])
SS_Cd_M_1000 = np.array([0.14,0.11,0.08,0.07,0.08,0.13,0.19,0.31])
SS_Cd_M_1500 = np.array([0.14,0.11,0.08,0.07,0.08,0.13,0.19,0.31])
SS_Cd_M_2000 = np.array([0.14,0.11,0.08,0.07,0.08,0.13,0.19,0.31])

def aerodynamics_coefficients(Mach,alpha):
    ''''
    The input Mach number is changed to the nearest Mach number for which aero data is available.
    Alpha is linearly interpolated between the datapoints, and linearly extrapolated afterwards.
    Function testing has been performed, for various Mach and alpha values.
    '''
    if Mach < 0.425:
        clvals = SS_Cl_M_0025
        cdvals = SS_Cd_M_0025
    if 0.425 <= Mach < 0.7:
        clvals = SS_Cl_M_0060
        cdvals = SS_Cd_M_0060
    if 0.7 <= Mach < 0.86:
        clvals = SS_Cl_M_0080
        cdvals = SS_Cd_M_0080
    if 0.86 <= Mach < 0.95:
        clvals = SS_Cl_M_0092
        cdvals = SS_Cd_M_0092
    if 0.95 <= Mach < 1.04:
        clvals = SS_Cl_M_0098
        cdvals = SS_Cd_M_0098
    if 1.04 <= Mach < 1.3:
        clvals = SS_Cl_M_0110
        cdvals = SS_Cd_M_0110
    if 1.3 <= Mach < 1.75:
        clvals = SS_Cl_M_0150
        cdvals = SS_Cd_M_0150
    if 1.75 <= Mach < 2.5:
        clvals = SS_Cl_M_0200
        cdvals = SS_Cd_M_0200
    if 2.5 <= Mach < 2.5:
        clvals = SS_Cl_M_0300
        cdvals = SS_Cd_M_0300
    if 3.5 <= Mach < 4.5:
        clvals = SS_Cl_M_0400
        cdvals = SS_Cd_M_0400
    if 4.5 <= Mach < 6.5:
        clvals = SS_Cl_M_0500
        cdvals = SS_Cd_M_0500
    if 6.5 <= Mach < 9:
        clvals = SS_Cl_M_0800
        cdvals = SS_Cd_M_0800
    if 9 <= Mach < 12.5:
        clvals = SS_Cl_M_1000
        cdvals = SS_Cd_M_1000
    if 12.5 <= Mach < 17.5:
        clvals = SS_Cl_M_1500
        cdvals = SS_Cd_M_1500
    if 17.5 <= Mach:
        clvals = SS_Cl_M_2000
        cdvals = SS_Cd_M_2000


    cd_function  = InterpolatedUnivariateSpline(SS_alpha_spacing,cdvals,k=1)
    cl_function = InterpolatedUnivariateSpline(SS_alpha_spacing,clvals,k=1)
    cl = cl_function(alpha)
    cd = cd_function(alpha)
    return cl,cd

