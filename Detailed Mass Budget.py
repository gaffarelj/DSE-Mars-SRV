# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:34:44 2020

@author: Dimitris
"""

from math import * 

k = "SSTO"

Isp       = 450
g         = 9.81
f         = 1
Ve        = Isp*g
TWR       = 1.5
Mass_caps = 6592.38
DV1       = 2844*f
DV2       = 3884*f 
 
def ClassIest(DeltaV,Ve,m_tot,TWR,m_upper):
    
    m_t     = Mass_caps    
    m_frac = exp(DeltaV/Ve)    
    m_prop = m_frac*m_tot/(1+m_frac)    
    m_RCS  = 0.0126*m_tot    
    Fvac   = TWR*3.7*m_tot    
    m_eng  = 0.00514*Fvac**0.92068    
    m_tank = 0.15*m_prop    
    m_thr_str = 1.949*10**(-3)*(Fvac/4.448)**1.0687*0.453
    m_stage = 0.001148*m_upper
    print(m_tank)

    rho_f = 423
    rho_ox = 1140
    Mf = m_prop/(1+1/3.8) #3.8 is F/o ratio
    Vf = Mf/rho_f
    Mox = m_prop - Mf
    Vox = Mox/rho_ox
    m_fuel_tank = (2.42-0.00271*rho_f*0.06242796047)*Vf*35.31466671**(0.8445+0.00047*rho_f*0.06242796047)*0.45359237
    m_oxid_tank = (2.42-0.00271*rho_ox*0.06242796047)*Vox*35.31466671**(0.8445+0.00047*rho_ox*0.06242796047)*0.45359237
    m_tank = m_oxid_tank+m_fuel_tank
    print(m_tank,m_prop)

    if m_stage == 0 : 
    
        M_total = m_t + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS
    
    if m_stage != 0 : 
    
        M_total =  m_tot + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS
        
    M_dry   = m_t + m_tank + m_eng +  m_thr_str
    
    return m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total, M_dry


def Class_I_spaceplane_est(DeltaV, Ve, m_tot, TWR, m_wing):
    m_t = Mass_caps
    m_frac = exp(DeltaV / Ve)
    m_prop = m_frac * m_tot / (1 + m_frac)
    m_RCS = 0.0126 * m_tot
    Fvac = TWR * 3.7 * m_tot
    m_eng = 0.00514 * Fvac ** 0.92068
    m_tank = 0.15 * m_prop
    m_thr_str = 1.949 * 10 ** (-3) * (Fvac / 4.448) ** 1.0687 * 0.453
    rho_f = 423
    rho_ox = 1140
    Mf = m_prop/(1+1/3.8) #3.8 is F/o ratio
    Vf = Mf/rho_f
    Mox = m_prop - Mf
    Vox = Mox/rho_ox
    m_fuel_tank = (2.42-0.00271*rho_f*0.06242796047)*Vf*35.31466671**(0.8445+0.00047*rho_f*0.06242796047)*0.45359237
    m_oxid_tank = (2.42-0.00271*rho_ox*0.06242796047)*Vox*35.31466671**(0.8445+0.00047*rho_ox*0.06242796047)*0.45359237
    m_tank = m_oxid_tank+m_fuel_tank
    M_total = m_t + m_prop + m_tank + m_thr_str + m_eng + m_RCS + m_wing

    M_dry = m_t + m_tank + m_eng + m_thr_str + m_wing

    return m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry,m_wing


if k == "SSTO" :
    m_upper = 0 
    DeltaV  = (DV1 + DV2)*f
    m_tot   = Mass_caps
    mt_margin = 0.1
    mts_margin = 0.1
    me_margin  = 0.1
    mc_margin    = 0.25 
    
    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total, M_dry  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper)
    
    print(m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total, M_dry)
    
    m_tot = M_total 
    
    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper)
    
    while M_total_new > M_total + 0.001 : 
        M_total = M_total_new 
        
        m_tot = M_total
        
        m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper)
        
        print(M_total_new,"--------------",M_total)
    
    #print(m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry)
    
    m_eng  = (1+me_margin) * m_eng 
    m_thr_str = m_thr_str*(1+mts_margin)
    m_tank = m_tank*(1+mt_margin) 
    m_caps = (mc_margin+1)* Mass_caps
    
    M_dry  = m_eng + m_thr_str + m_tank + m_caps
    
if k =="2_stage" : 
    m_upper = 0 
    m_tot   = Mass_caps
    mt_margin = 0.1
    mts_margin = 0.1
    me_margin  = 0.1
    mc_margin    = 0.25 
    
    # 1st stage calculation 
    
    m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1, M_dry1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper)
    
    m_upper = M_total1 
    m_tot   = M_total1 
    
    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total, M_dry  = ClassIest(DV2,Ve,m_tot,TWR,m_upper)
    
    m_tot = M_total1
    m_upper = 0 
    
    m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1, M_dry1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper)
    
    m_upper = M_total1 
    m_tot   = M_total1 
    
    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry  = ClassIest(DV2,Ve,m_tot,TWR,m_upper)
    
    
    while M_total_new > M_total + 0.001 : 
        M_total = M_total_new 
        m_tot   = M_total1 
        m_upper = 0 
        
        m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1, M_dry1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper)
    
        m_upper = M_total1 
        m_tot   = M_total1 
        
        m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry  = ClassIest(DV2,Ve,m_tot,TWR,m_upper)
        
        print(M_total_new,"-------",M_total1)

def takeoff_wing_sizing_shuttle_like(M_takeoff,takeoff_mach,takeoff_cl):
    rho = 0.02 #Density at surface on Mars
    vel = 240 * takeoff_mach #takeoff speed as function of speed of sound (240) in m/s
    S = M_takeoff*3.7/(0.5*rho*vel**2*takeoff_cl)
    b = 30.5/874.5* 8 * S #Adjusted manually from spaceshuttle planform
    taper = 0.2 #Roughly, from spaceshuttle planform again
    c_root = S / (0.5*b)*(1/(1+taper))
    c_tip = taper * c_root
    return S, b, c_root, c_tip


def AVID_wing_mass( M_land , b , S_exp , c_root , tc ):
    M_land = M_land * 2.205
    b = b / 0.3048
    print("Wing Area:", S_exp)
    S_exp = S_exp / (0.3048**2) * 0.8 #90% of wing exposed

    c_root = c_root / 0.3048

    #Mwing = 1575*((M_land*3.75*b*S_exp)/(c_root*tc*10**9))**0.67
    #Mwing = Mwing * 0.4536

    Mwing = 1.498 * S_exp **1.176 * 0.4536 * 0.2
    print("Wing Mass:", Mwing)

    return Mwing

if k == "SPACEPLANE":
    m_upper = 0
    TWR = 0.15
    DeltaV = (DV1 + DV2) * f - 1000 # Assumed no DeltaV for landing
    m_tot = Mass_caps
    mt_margin = 0.1
    mts_margin = 0.1
    me_margin = 0.1
    mc_margin = 0.25
    m_wing = 1000

    m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry, m_wing = Class_I_spaceplane_est(DeltaV,Ve,m_tot,TWR,m_wing)

    print(m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry, m_wing)

    S,b,c_root,c_tip = takeoff_wing_sizing_shuttle_like(M_total,0.8,1.5)

    m_wing = AVID_wing_mass(M_dry,b,S,c_root,0.15)

    m_tot = M_total

    m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total_new, M_dry, m_wing = Class_I_spaceplane_est(
        DeltaV, Ve, m_tot, TWR, m_wing)

    while M_total_new > M_total + 0.001:
        M_total = M_total_new

        S, b, c_root, c_tip = takeoff_wing_sizing_shuttle_like(M_total, 0.8, 1.5)
        m_wing = AVID_wing_mass(M_dry,b,S, c_root, 0.15)

        m_tot = M_total

        m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total_new, M_dry, m_wing = Class_I_spaceplane_est(
            DeltaV, Ve, m_tot, TWR, m_wing)
        print(m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry, m_wing)
        #print(M_total_new,m_wing,"--------------",M_total)

    # print(m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry)

    m_eng = (1 + me_margin) * m_eng
    m_thr_str = m_thr_str * (1 + mts_margin)
    m_tank = m_tank * (1 + mt_margin)
    m_caps = (mc_margin + 1) * Mass_caps

    M_dry = m_eng + m_thr_str + m_tank + m_caps


