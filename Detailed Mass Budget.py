# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:34:44 2020

@author: Dimitris
"""

from math import * 

k = "SSTO"

N_crew = 6
N_days = 4
conv   = 0.4535
fcom   = 0.25
fEPs   = 0.25
fls    = 0.4
Klg    = 0.033

m_fuel_cells   = (3030*N_crew/7)*conv
m_bat          = 216

m_EPS          = m_fuel_cells + m_bat

m_Life_support = (2444*N_crew/7 + 645*N_crew + 86.4*N_days)*conv

m_Com          = (131*N_days/7 + 1400*N_crew/7)*conv
m_data_hand    = (302 + 828*N_days/7 + 1010*N_crew/7)*conv 
m_GNC          = (242 + 108*N_days/7 + 617*N_crew/7)*conv 

m_av           = m_GNC + m_data_hand + m_Com*(1-fcom)



m_cabin       = (28.31*(39.66*(N_crew*N_days)**1.002)**0.6916)*conv
#m_pl_doors    = (0.257*barea/2)* conv


m_pl           = 1000
m_payl_cont    = 700


Isp       = 410 
g         = 9.81
f         = 1
Ve        = Isp*g
TWR       = 1.5
Mass_caps = m_EPS*(1-fEPs) + m_Life_support*(1-fls) + m_av + m_cabin + m_pl + m_payl_cont
DV1       = 2200*f
DV2       = 3200*f

mlg_margin  = 0.2
mt_margin = 0.2
mts_margin = 0.2
me_margin  = 0.2
mc_margin    = 0.2

def ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2):
    
    klg     = 0.033
    m_t     = Mass_caps
    m_frac = exp(DeltaV/Ve)    
    m_prop = m_frac*m_tot/(1+m_frac)
    m_stage = 0.001148*m_upper


    rho_f = 423
    rho_ox = 1140
    Mf = m_prop/(1+1/3.8) #3.8 is F/o ratio
    Vf = Mf/rho_f
    Mox = m_prop - Mf
    Vox = Mox/rho_ox
    #m_fuel_tank = (2.42-0.00271*rho_f*0.06242796047)*(Vf*35.31466671)**(0.8445+0.00047*rho_f*0.06242796047)*0.45359237
    #m_oxid_tank = (2.42-0.00271*rho_ox*0.06242796047)*(Vox*35.31466671)**(0.8445+0.00047*rho_ox*0.06242796047)*0.45359237
    #m_tank = m_oxid_tank+m_fuel_tank
    #print(m_tank,m_fuel_tank,m_oxid_tank,m_prop)
    #print(m_tank)

    m_tank = (Vf+Vox)*3*10**6/(6.43*10**4) #3MPa assumed MEOP, based on a reference

    if m_stage == 0 : 
    
        Fvac   = TWR*3.7*m_tot
        m_RCS  = 0.0126*m_tot
        m_eng  = 0.00514*Fvac**0.92068
        #m_eng  = 1.866**(-10)*Fvac + 0.00130*Fvac+77.4
        m_thr_str = (1.949*10**(-3)*(Fvac/4.448)**1.0687)*0.453
        M_dry     = m_eng + m_thr_str + m_tank + Mass_caps
        Mlg       = klg * M_dry*1.15
        M_total = m_t + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS + Mlg
    
    if m_stage != 0 :
        M_dry   = m_t + m_tank + m_eng +  m_thr_str

    return m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total, M_dry


def Class_I_spaceplane_est(DeltaV, Ve, m_tot, TWR, m_wing):
    m_t = Mass_caps
    m_frac = exp(DeltaV / Ve)
    m_prop = m_frac * m_tot / (1 + m_frac)
    m_RCS = 0.0126 * m_tot
    Fvac = TWR * 3.7 * m_tot
    m_eng = 0.00514 * Fvac ** 0.92068
    m_thr_str = 1.949 * 10 ** (-3) * (Fvac / 4.448) ** 1.0687 * 0.453
    m_landinggear = 0.010784* ((m_tot-m_prop-m_RCS)*0.453)**1.0861 * 0.453
    rho_f = 423
    rho_ox = 1140
    Mf = m_prop/(1+1/3.8) #3.8 is F/o ratio
    Vf = Mf/rho_f
    Mox = m_prop - Mf
    Vox = Mox/rho_ox
    #m_fuel_tank = (2.42-0.00271*rho_f*0.06242796047)*(Vf*35.31466671)**(0.8445+0.00047*rho_f*0.06242796047)*0.45359237
    #m_oxid_tank = (2.42-0.00271*rho_ox*0.06242796047)*(Vox*35.31466671)**(0.8445+0.00047*rho_ox*0.06242796047)*0.45359237
    #m_tank = m_oxid_tank+m_fuel_tank
    Fvac   = TWR*3.7*m_tot
    m_RCS  = 0.0126*m_tot
    m_eng  = 0.00514*Fvac**0.92068
    m_thr_str = 1.949*10**(-3)*(Fvac/4.448)**1.0687*0.453
    #M_total =  m_tot + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS

    #return m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total
    #M_dry     = m_eng + m_thr_str + m_tank + Mass_caps
    #Mlg       = klg * M_dry*1.15
    #M_total =  m_tot + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS + Mlg

    #return m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total,Mlg

    #print(m_tank)

    m_tank = (Vf+Vox)*3*10**6/(6.43*10**4)
    #print(m_tank, m_prop)
    #print("-----")
    M_total = m_t + m_prop + m_tank + m_thr_str + m_eng + m_RCS + m_wing + m_landinggear

    M_dry = m_t + m_tank + m_eng + m_thr_str + m_wing + m_landinggear

    return m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry,m_wing


if k == "SSTO" :
    m_tot2 =  0 
    m_upper = 0 
    DeltaV  = (DV1 + DV2)*f
    DeltaV  = 5400
    m_tot   = Mass_caps
    
    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total,Mlg  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2)
    
    #print(m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total)
    
    m_tot = M_total 
    
    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2)
    
    while M_total_new > M_total + 0.001 : 
        M_total = M_total_new 
        
        m_tot = M_total
        
        m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2)
        
        print(M_total_new,"--------------",M_total)
    
    #print(m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new, M_dry)
    
    m_eng  = (1+me_margin) * m_eng 
    m_thr_str = m_thr_str*(1+mts_margin)
    m_tank = m_tank*(1+mt_margin) 
    m_caps = (mc_margin+1)* Mass_caps
    Mlg    = (1+mlg_margin)*Mlg

    M_dry  = m_eng + m_thr_str + m_tank + m_caps + Mlg
    
    print(M_dry,m_prop,M_total,m_tank,m_thr_str,m_eng)
    
    
if k =="2_stage" : 
    m_upper = 0 
    m_tot   = Mass_caps
    m_tot2    = 0
    
    # 1st stage calculation 
    
    m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1,Mlg1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper,m_tot2)
    
    m_upper = M_total1 
    m_tot   = M_total1 
    m_tot2    = 0 

    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total,Mlg  = ClassIest(DV2,Ve,m_tot,TWR,m_upper,m_tot2)
    
    m_tot = M_total1
    m_upper = 0 
    m_tot2    = M_total 

    m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1,Mlg1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper,m_tot2)
    
    m_upper = M_total1 
    m_tot   = M_total1 

    m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg  = ClassIest(DV2,Ve,m_tot,TWR,m_upper,m_tot2)

    
    
    while M_total_new > M_total + 0.001 : 
        M_total = M_total_new 
        m_tot   = M_total1 
        m_tot2    = M_total_new 
        m_upper = 0 
        
        m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1,Mlg1= ClassIest(DV1,Ve,m_tot,TWR,m_upper,m_tot2)
    
        m_upper = M_total1 
        m_tot   = M_total1 
        
        m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg  = ClassIest(DV2,Ve,m_tot,TWR,m_upper,m_tot2)
        
        print(M_total_new,"-------",M_total)
        
    m_eng  = (1+me_margin) * m_eng 
    m_thr_str = m_thr_str*(1+mts_margin)
    m_tank = m_tank*(1+mt_margin) 
    m_caps = (mc_margin+1)* Mass_caps
    Mlg    = (1+mlg_margin)*Mlg

    M_dry  = m_eng + m_thr_str + m_tank 
    
    m_eng1  = (1+me_margin) * m_eng1 
    m_thr_str1 = m_thr_str1*(1+mts_margin)
    m_tank1 = m_tank1*(1+mt_margin) 
    Mlg1    = (1+mlg_margin)*Mlg1


    M_dry1  = m_eng1 + m_thr_str1 + m_tank1 + Mlg1

    M_dry_tot = M_dry1 + M_dry + m_caps + Mlg
    
    print("---------------------2st stage-------------")
    print("Tank stage 2", m_tank1)
    print("Thrust Stuct stage 2", m_thr_str1)
    print("eng stage 2", m_eng1)
    print("Dry stage 2", M_dry1)
    print("---------------------1st stage-------------")
    print("Tank stage 1", m_tank)
    print("Thrust Stuct stage 1", m_thr_str)
    print("Eng stage 1", m_eng)
    print("Dry stage 1", M_dry)
    print("Capsule Mass", m_caps)
    print("Total Mass (Dry):",M_dry_tot)
    print("Total Mass ",M_total)
    print("Total propellant mass",m_prop1 + m_prop)
    print("propellant stage 2",m_prop1)
    print("propellant stage 1",m_prop)

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

    Mwing = 1.498 * S_exp **1.176 * 0.4536 * 0.4
    print("Wing Mass:", Mwing)

    return Mwing

if k == "SPACEPLANE":
    m_upper = 0
    TWR = 1.5
    DeltaV = 5400 * f - 500 # Assumed no DeltaV for landing
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

        S, b, c_root, c_tip = takeoff_wing_sizing_shuttle_like(M_dry*1.1, 0.8, 1.5)
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


