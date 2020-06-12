import numpy as np
from matplotlib import pyplot as plt

def gravitygradient_disturbance(Iyy,Izz,omega,theta):
    #theta: max deviation between nadir and local vertical on XZ-plane
    #phi: max deviation between nadir and local vertical on XY-plane
    #n: orbital mean motion
    # n = np.sqrt(mu/(2*a**3))

    Tg = 3*omega**2*(Izz-Iyy)*np.sin(2*theta)
    # Tg = np.array((Tg_x,Tg_y,0))
    return Tg

def Drag_surface_cp(alpha):
    length_engine  = 3.30 #m
    width_engine   = 1.43 #m
    chord_wings   = 3     #m
    span_wings    = 5     #m
    length_body    = 14.01#m
    width_body     = 7.61 #m
    radius_nose    = 3.6  #m

    A_engines      = np.sin(alpha) * 3 * 1.43/2*3.30                  #m^2
    A_wings        = np.sin(alpha) * 2 * chord_wings * span_wings     #m^2
    A_body         = np.sin(alpha) * 14.01*7.610                      #m^2
    A_nose         = np.cos(alpha) * np.pi* radius_nose**2            #m^2
    A_tot          = A_nose + A_engines + A_body                      #m^2

    d_engines      = length_engine / 2                                                     #m
    d_wings        = length_engine + chord_wings / 2                                       #m
    d_body         = length_engine + length_body / 2                                       #m
    d_nose         = length_engine + length_body + 4 * radius_nose / (3 * np.pi)           #m

    Cp             = (A_engines * d_engines + A_body * d_body + A_nose * d_nose) / (A_tot) #m

    return A_tot, Cp


def Drag_force(q,CD,S):
    drag = q*CD*S
    return drag

def aerodynamic_disturbance(Cp,Cg,Drag,alpha):
    #Cg: location of center of gravity
    #Cp: location of center of pressure
    r = (Cp-Cg) * np.sin(alpha)
    Td = Drag * r
    return Td

def solarpressure_disturbance(alpha,Cg):
    length_engine  = 3.30 #m
    width_engine   = 1.43 #m
    chord_wings    = 5    #m
    span_wings     = 3    #m
    length_body    = 14.01#m
    width_body     = 7.61 #m
    radius_nose    = 3.6  #m

    A_engines      = np.cos(alpha) * 3 * 1.43/2*3.30                  #m^2
    A_wings        = np.cos(alpha) * 2 * chord_wings * span_wings     #m^2
    A_body         = np.cos(alpha) * 14.01*7.610                      #m^2
    A_nose         = np.cos(alpha) * np.pi* radius_nose**2            #m^2
    A_tot          = A_nose + A_engines + A_body + A_wings            #m^2

    d_engines      = length_engine / 2                                                     #m
    d_wings        = length_engine + chord_wings / 2                                       #m
    d_body         = length_engine + length_body / 2                                       #m
    d_nose         = length_engine + length_body + 4 * radius_nose / (3 * np.pi)           #m
    Cp             = (A_engines * d_engines + A_body * d_body + A_nose * d_nose) / (A_tot) #m

    As = A_tot
    #Fs: Solar irradance on Mars surface [W/m^2]
    #c:  Speed of light
    Fs = 586.2
    c  = 3. * 10**8
    q  = 0.6
    F = Fs/c*As*(1+q)*np.cos(alpha)
    Tsp = F*(Cp-Cg)
    return Tsp

def magnetic_disturbance(R):
    #M: Magnetic moment of Mars [tesla*m^3]
    #B: Magnetic field of Mars at orbit height
    #D: Residual dipole of vehicle [A*m^2]
    D = 0.1
    M = 1.22*10**12
    B = 2*M/R**3
    Tm = D*B
    return Tm
print(magnetic_disturbance(3389.5*10**3))
print(gravitygradient_disturbance)
