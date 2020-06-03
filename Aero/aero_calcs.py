#Made by the aerodynamicist

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import interpolate

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

cd_data_SS = np.vstack([SS_Cd_M_0025,SS_Cd_M_0060,SS_Cd_M_0080,SS_Cd_M_0092,SS_Cd_M_0098,SS_Cd_M_0110,SS_Cd_M_0150,SS_Cd_M_0200,SS_Cd_M_0300,SS_Cd_M_0400,SS_Cd_M_0500,SS_Cd_M_0800,SS_Cd_M_1000,SS_Cd_M_1500,SS_Cd_M_2000])
cl_data_SS = np.vstack([SS_Cl_M_0025,SS_Cl_M_0060,SS_Cl_M_0080,SS_Cl_M_0092,SS_Cl_M_0098,SS_Cl_M_0110,SS_Cl_M_0150,SS_Cl_M_0200,SS_Cl_M_0300,SS_Cl_M_0400,SS_Cl_M_0500,SS_Cl_M_0800,SS_Cl_M_1000,SS_Cl_M_1500,SS_Cl_M_2000])
f_cl_SS = interpolate.interp2d(SS_alpha_spacing, [0.25, 0.6, 0.8, 0.92, 0.98, 1.1, 1.5, 2, 3, 4, 5, 8, 10, 15, 20], cl_data_SS, kind='cubic')
f_cd_SS = interpolate.interp2d(SS_alpha_spacing, [0.25, 0.6, 0.8, 0.92, 0.98, 1.1, 1.5, 2, 3, 4, 5, 8, 10, 15, 20], cd_data_SS, kind='cubic')


def SS_aerodynamics_coefficients(mach,alpha):
    cl = f_cl_SS(alpha, mach)
    cd = f_cd_SS(alpha, mach)

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

cd_data_H = np.vstack([H_Cd_M_0120,H_Cd_M_0150,H_Cd_M_0200,H_Cd_M_0300,H_Cd_M_0500,H_Cd_M_1000,H_Cd_M_2000])
cl_data_H = np.vstack([H_Cl_M_0120,H_Cl_M_0150,H_Cl_M_0200,H_Cl_M_0300,H_Cl_M_0500,H_Cl_M_1000,H_Cl_M_2000])
f_cl_H = interpolate.interp2d(H_alpha_spacing, [1.2, 1.5, 2, 3, 5, 10, 20], cl_data_H, kind='cubic')
f_cd_H = interpolate.interp2d(H_alpha_spacing, [1.2, 1.5, 2, 3, 5, 10, 20], cd_data_H, kind='cubic')

def H_aerodynamics_coefficients(mach,alpha):
    cl = f_cl_H(alpha, mach)/2
    cd = f_cd_H(alpha, mach)

    if alpha > 45:
        cd = 1.22
        cl = 0.2

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
