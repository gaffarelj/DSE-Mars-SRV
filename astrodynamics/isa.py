from math import *

def ISA(h):
    To = 288.15
    g = 9.81
    R = 287
    po = 101325
    a1  =  -0.0065    #K/m                                                                  #Both type layers, Density: rh= p/(RT)
    a2  =   0         #K/m
    a3  =   0.001     #K/m
    a4  =   0.0028    #K/m
    a5  =   0         #K/m
    a6  =  -0.0028    #K/m
    a7  =  -0.002     #K/m

    if h>0:
        To = To #K
        po=101325   #Pa
        h1= min(h,11000)                    #Troposhpere 1, Gradient layer
        T=To+h1*a1
        p=po*(T/To)**(-g/(R*a1))
        rh=p/(R*T)
    if h>11000:
        To=T #K
        po=p #Pa
        h2=min(h,20000)
        T=To+(h2-h1)*a2                     #Troposphere 2, Isothermal Layer
        p=po*exp(-g*(h2-h1)/(R*T))
        rh=p/(R*T)
    if h>20000:
        To=T #K
        po=p #Pa
        h3=min(h,32000)
        T=To+(h3-h2)*a3                     #Stratosphere 1, Gradient Layer
        p=po*(T/To)**(-g/(R*a3))
        rh=p/(R*T)
    if h>32000:
        To=T #K                                                                                                   #ISA CALCULATOR, All layers
        po=p #Pa
        h4=min(h,47000)
        T=To+(h4-h3)*a4                     #Stratosphere 2, Gradient Layer
        p=po*(T/To)**(-g/(R*a4))
        rh=p/(R*T)
    if h>47000:
        To=T #K
        po=p #Pa
        h5=min(h,51000)
        T=To+(h5-h4)*a5                     #Stratosphere 2, Isothermal Layer
        p=po*exp(-g*(h5-h4)/(R*T))
        rh=p/(R*T)
    if h>51000:
        To=T #K
        po=p #Pa
        h6=min(h,71000)
        T=To+(h6-h5)*a6                     #Mesosphere 1, Gradient Layer
        p=po*(T/To)**(-g/(R*a6))
        rh=p/(R*T)
    if h>71000:
        To=T #K
        po=p #Pa
        h7=min(h,86000)                     #Mesosphere 1, Gradient Layer
        T=To+(h7-h6)*a7
        p=po*(T/To)**(-g/(R*a7))
        rh=p/(R*T)

    return ( T, p , rh ) 

