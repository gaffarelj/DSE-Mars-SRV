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


def SS_aerodynamics_coefficients(Mach,alpha):
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
        clvals = SS_Cl_M_0200
        cdvals = SS_Cd_M_0200
    if 2.5 <= Mach < 3.5:
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


H_alpha_spacing = np.array([0,5,10,15,20,25,30,35,40,45])

H_Cl_M_0120 = np.array([-0.02,0.27,0.56,0.86,1.17,1.17,1.17,1.17,1.17,1.17])
H_Cl_M_0150 = np.array([-0.02,0.16,0.36,0.57,0.79,1.01,1.21,1.38,1.38,1.38])
H_Cl_M_0200 = np.array([-0.03,0.10,0.26,0.42,0.60,0.78,0.94,1.08,1.20,1.28])
H_Cl_M_0300 = np.array([-0.04,0.06,0.19,0.33,0.48,0.62,0.77,0.90,1.01,1.08])
H_Cl_M_0500 = np.array([-0.04,0.05,0.14,0.26,0.40,0.54,0.68,0.80,0.91,0.99])
H_Cl_M_1000 = np.array([-0.05,0.03,0.12,0.23,0.36,0.50,0.63,0.75,0.86,0.95])
H_Cl_M_2000 = np.array([-0.05,0.01,0.07,0.16,0.28,0.40,0.53,0.66,0.77,0.85])

H_Cd_M_0120 = np.array([0.10,0.12,0.20,0.33,0.52,0.77,1.08,1.08,1.08,1.08])
H_Cd_M_0150 = np.array([0.09,0.10,0.14,0.23,0.36,0.54,0.77,1.03,1.32,1.32])
H_Cd_M_0200 = np.array([0.08,0.09,0.12,0.18,0.29,0.42,0.60,0.82,1.07,1.33])
H_Cd_M_0300 = np.array([0.07,0.08,0.10,0.15,0.23,0.35,0.51,0.69,0.90,1.21])
H_Cd_M_0500 = np.array([0.07,0.07,0.09,0.14,0.20,0.31,0.45,0.63,0.82,1.05])
H_Cd_M_1000 = np.array([0.07,0.07,0.08,0.12,0.18,0.29,0.42,0.60,0.79,1.02])
H_Cd_M_2000 = np.array([0.05,0.05,0.06,0.09,0.14,0.24,0.36,0.52,0.70,0.92])

def H_aerodynamics_coefficients(Mach,alpha):
    ''''
plt.plot(SS_alpha_spacing,SS_Cd_M_0500)
plt.show()
plt.plot(SS_alpha_spacing,SS_Cl_M_0500)
plt.show()
e datapoints, and linearly extrapolated afterwards.
    Function testing has been performed, for various Mach and alpha values.
    '''
    if Mach < 1.25:
        clvals = H_Cl_M_0120
        cdvals = H_Cd_M_0120
    if 1.25 <= Mach < 1.75:
        clvals = H_Cl_M_0150
        cdvals = H_Cd_M_0150
    if 1.75 <= Mach < 2.5:
        clvals = H_Cl_M_0200
        cdvals = H_Cd_M_0200
    if 2.5 <= Mach < 4:
        clvals = H_Cl_M_0300
        cdvals = H_Cd_M_0300
    if 4 <= Mach < 7.5:
        clvals = H_Cl_M_0500
        cdvals = H_Cd_M_0500
    if 7.5 <= Mach < 15:
        clvals = H_Cl_M_1000
        cdvals = H_Cd_M_1000
    if 15 <= Mach:
        clvals = H_Cl_M_2000
        cdvals = H_Cd_M_2000

    cd_function  = InterpolatedUnivariateSpline(H_alpha_spacing,cdvals,k=3)
    cl_function = InterpolatedUnivariateSpline(H_alpha_spacing,clvals,k=3)
    cl = cl_function(alpha)/2
    cd = cd_function(alpha)
    return cl,cd


if __name__ == '__main__':
    plt.plot(SS_alpha_spacing,SS_Cd_M_0500)
    plt.show()
    plt.plot(SS_alpha_spacing,SS_Cl_M_0500)
    plt.show()

    plt.plot(H_alpha_spacing,H_Cd_M_0500)
    plt.show()
    plt.plot(H_alpha_spacing,H_Cl_M_0500)
    plt.show()
