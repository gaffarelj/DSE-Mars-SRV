# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:34:44 2020

@author: Dimitris
"""

from math import * 

k = "2_stage"

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
    
    if m_stage == 0 : 
    
        M_total = m_t + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS
    
    if m_stage != 0 : 
    
        M_total =  m_tot + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS
        
    M_dry   = m_t + m_tank + m_eng +  m_thr_str
    
    return m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total, M_dry

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
        
        #print(M_total_new,"--------------",M_total)
    
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
        

    