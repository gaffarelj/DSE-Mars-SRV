# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:06:20 2020

@author: lucat
"""

import numpy as np

#========================================================================================================================================================================================================================
# Functions 
#========================================================================================================================================================================================================================
#Mach
def M(V,a):
    M=V/a
    return M
#speed of sound
def a(gamma,R,T):
    a=np.sqrt(gamma*Rgas*T)
    return a
#mu=gravitational parameter Mars, Req=equatorial radius Mars, R_mars=mean volumetric radius Mars, h=altitude, theta=latitude [rad]
#gravitational acceleration
def g(mu,Req,R,h,theta):
    P2=3/2*np.sin(theta)**2-1/2
    g=mu/(R+h)**2*(1-3*J2*(Req/(R+h))**2*P2)
    return g

#From mars_standard_atmosphere code [based on Viking 1 measurements]:

#pressure
def p(h):

    if h < 39000:
        p = 610.5*(228.5/(228.5-1.8*h/1000))**(19.435/-1.8)
    elif h < 48000:
        p = 11.6025*np.exp(-19.435*(h/1000-39)/158.3)
    elif h < 55000:
        p = 3.84305*(158.3/(158.3-2.35*(h/1000-48)))**(19.435/-2.35)
    elif h < 66000:
        p = 1.55091*(141.85/(141.85+0.65*(h/1000-55)))**(19.435/0.65)
    elif h < 75000:
        p = 0.356464*(149/(149-2.5*(h/1000-66)))**(19.435/-2.5)
    elif h < 84000:
        p = 0.099843*(126.5/(126.5+2.5*(h/1000-75)))**(19.435/2.5)
    elif h < 95000:
        p = 0.0279653*np.exp(-19.435*(h/1000-84)/149)
    elif h < 105000:
        p = 0.00666032*(149/(149-1.4*(h/1000-95)))**(19.435/-1.4)
    elif h < 115897:
        p = 0.00169282*(135/(135-0.65*(h/1000-105)))**(19.435/-0.65)
    else:
        p = 0
    return p 

#Temperature
def T(h):

    if h < 39000:
        T = 228.5 - 1.8*h/1000
    elif h < 48000:
        T = 158.3
    elif h < 55000:
        T = 271.1 - 2.35*h/1000
    elif h < 66000:
        T = 106.1 + 0.65*h/1000
    elif h < 75000:
        T = 314 - 2.5*h/1000
    elif h < 84000:
        T = -61 + 2.5*h/1000
    elif h < 95000:
        T = 149
    elif h < 105000:
        T = 282 - 1.4*h/1000
    elif h < 115897:
        T = 203.25 - 0.65*h/1000
    else:
        e = (h/1000 - 120)*(3389.51 + 120)/(3389.51 + h/1000)
        T = 200 - 72.225*np.exp(-0.0195*e)
    return T

#density
def rho(p, T, Rgas):
    return p/(Rgas*T)

#equivalent jet velocity
def ceff(Isp):
    ceff=9.80665*Isp
    return ceff

#Mass ratio. Note M0=Mi+Mp
def Lambda(M0,Me):
    Lambda=M0/Me
    return Lambda


#propellant mass: rewritten Tsiolkovksy in terms of wet mass
def Mprop(ceff,Mwet,deltaV):
    Mprop=(Mwet*(np.exp(deltaV/ceff)-1))/np.exp(deltaV/ceff)
    return Mprop




