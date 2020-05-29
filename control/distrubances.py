import numpy as np
from matplotlib import pyplot as plt 

def gravitygradient_disturbance(Ixx,Iyy,Izz,a,phi,theta):
    #theta: max deviation between nadir and local vertical on XZ-plane
    #phi: max deviation between nadir and local vertical on XY-plane
    #n: orbital mean motion
    n = np.sqrt(mu/(2*a**3))
    Tg_x = 3*n**2*(Izz-Iyy)*np.sin(2*phi)
    Tg_y = 3*n**2*(Izz-Ixx)*np.sin(2*theta)
    Tg = np.array((Tg_x,Tg_y,0))
    return Tg

def aerodynamic_disturbance(Cg,Cp,Drag):
    #Cm: location of center of mass
    #Cp: location of center of pressure
    #Drag = 0.5*rho*V**2*S*CD
    r = Cg - Cp
    Ta = np.cross(Drag,r)
    return Ta

def solarpressure_disturbance(Fs):
    #Fs: Solar irradance on Mars surface [W/m^2]
    Fs = 586.2
    F = Fs/c*As*(1+q)*cos(I_angle)
    Tsp = F*(Cps-Cg)

def magnetic_disturbance():
    #M: Magnetic moment of Mars [tesla*m^3]
    #B: Magnetic field of Mars at orbit height
    #D: Residual dipole of vehicle [A*m^2]
    D = 0.1
    M = 1.22*10**12
    B = 2*M/R**3
    Tm = D*B