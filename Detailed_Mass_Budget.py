# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:34:44 2020

@author: Dimitris
"""

from math import *
import numpy as np

def Mass_conc(DV1,DV2,Isp,k,MGA="YES") :
	N_crew = 6
	N_days = 4
	conv   = 0.4535
	
	if MGA == "YES":
	
		MGAEPS          = 0.117
		MGA_data_hand   = 0.25 
		MGATPS          = 0.25
		MGA_com         = 0.461
		MGA_GNC         = 0.215
		MGA_cab         = 0.375
		MGA_ls          = 0.225
		MGA_se_power    = 0.5
	
	else :
	
		MGAEPS  = 0
		MGA_data_hand   = 0
		MGATPS  = 0
		MGA_com = 0
		MGA_GNC = 0
		MGA_cab = 0
		MGA_ls  = 0
		MGA_se_power    = 0
	
	m_fuel_cells   = (3030*N_crew/7)*conv
	m_bat          = 216
	
	m_EPS          = m_fuel_cells + m_bat
	
	m_Life_support = (2444*N_crew/7 + 645*N_crew + 86.4*N_days)*conv
	
	m_Com          = (131*N_days/7 + 1400*N_crew/7)*conv
	m_data_hand    = (302 + 828*N_days/7 + 1010*N_crew/7)*conv 
	m_GNC          = (242 + 108*N_days/7 + 617*N_crew/7)*conv 
	
	m_av           = m_GNC*(1+MGA_GNC) + m_data_hand*(1+MGA_data_hand) + m_Com*(1+MGA_com)
	
	M_TPS          = 2500
	
	m_cabin        = (28.31*(39.66*(N_crew*N_days)**1.002)**0.6916)*conv
	
	m_pl           = 1200
	m_payl_cont    = 0.7*m_pl
	
	
	Isp       = Isp  
	g         = 9.81
	f         = 1
	Ve        = Isp*g
	TWR       = 1.5
	Mass_caps = m_EPS*(1+MGAEPS) + m_Life_support*(1+MGA_ls) + m_av + m_cabin*(1+MGA_cab) + m_pl + m_payl_cont + M_TPS*(1+MGATPS)
	DV1       = DV1*f
	DV2       = DV2*f
	
	def ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2):
		
		klg         = 0.033
		m_t         = Mass_caps
		m_frac      = exp(DeltaV/Ve)    
		m_prop      = m_frac*m_tot/(1+m_frac)
		m_stage     = 0.001148*m_upper
	
	
		rho_f   = 423
		rho_ox  = 1140
		Mf      = m_prop/(1+1/3.8) #3.8 is F/o ratio
		Vf      = Mf/rho_f
		Mox     = m_prop - Mf
		Vox     = Mox/rho_ox

		m_tank  = (Vf+Vox)*3*10**6/(6.43*10**4) #3MPa assumed MEOP, based on a reference
	
		if m_stage == 0 : 
		
			Fvac        = TWR*3.7*m_tot
			m_RCS       = 0.0126*m_tot
			m_eng       = 0.00514*Fvac**0.92068
			m_thr_str   = (1.949*10**(-3)*(Fvac/4.448)**1.0687)*0.453
			M_dry       = m_eng + m_thr_str + m_tank + Mass_caps
			Mlg         = klg * M_dry*1.15
			M_total     = m_t + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS + Mlg
		
		if m_stage != 0 :
	
			Fvac        = TWR*3.7*m_tot2
			m_RCS       = 0.0126*m_tot2
			m_eng       = 0.00514*Fvac**0.92068
			m_thr_str   = 1.949*10**(-3)*(Fvac/4.448)**1.0687*0.453
			M_dry       = m_eng + m_thr_str + m_tank + Mass_caps
			Mlg         = klg * M_dry*1.15
			M_total     =  m_tot + m_prop + m_tank + m_thr_str + m_eng + m_stage + m_RCS + Mlg
	
		return m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total,Mlg
	
	def Class_I_spaceplane_est(DeltaV, Ve, m_tot, TWR, m_wing):
		
		m_t             = Mass_caps
		m_frac          = exp(DeltaV / Ve)
		m_prop          = m_frac * m_tot / (1 + m_frac)
		m_RCS           = 0.0126 * m_tot
		Fvac            = TWR * 3.7 * m_tot
		m_eng           = 0.00514 * Fvac ** 0.92068
		m_thr_str       = 1.949 * 10 ** (-3) * (Fvac / 4.448) ** 1.0687 * 0.453
		m_landinggear   = 0.010784* ((m_tot-m_prop-m_RCS)*0.453)**1.0861 * 0.453
		rho_f           = 423
		rho_ox          = 1140
		Mf              = m_prop/(1+1/3.8) #3.8 is F/o ratio
		Vf              = Mf/rho_f
		Mox             = m_prop - Mf
		Vox             = Mox/rho_ox

		Fvac      = TWR*3.7*m_tot
		m_RCS     = 0.0126*m_tot
		m_eng     = 0.00514*Fvac**0.92068
		m_thr_str = 1.949*10**(-3)*(Fvac/4.448)**1.0687*0.453
	
		m_tank    = (Vf+Vox)*3*10**6/(6.43*10**4)

		M_total   = m_t + m_prop + m_tank + m_thr_str + m_eng + m_RCS + m_wing + m_landinggear
	
		M_dry     = m_t + m_tank + m_eng + m_thr_str + m_wing + m_landinggear

		return m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry,m_wing
	

	if k == "SSTO" :
		m_tot2  =  0 
		m_upper = 0 
		DeltaV  = (DV1 + DV2)*f
		m_tot   = Mass_caps
		
		m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total,Mlg  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2)
			   
		m_tot = M_total 
		
		m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2)
		
		while M_total_new > M_total + 0.001 : 
			M_total = M_total_new 
			
			m_tot   = M_total
			
			m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg  = ClassIest(DeltaV,Ve,m_tot,TWR,m_upper,m_tot2)
			 
		M_dry  = m_eng + m_thr_str + m_tank + Mass_caps + Mlg
		
		return [M_dry,m_prop]
		
		
	if k =="2_stage" : 
		m_upper     = 0 
		m_tot       = Mass_caps
		m_tot2      = 0
		
		# 1st stage calculation 
		
		m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1,Mlg1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper,m_tot2)
		
		m_upper     = M_total1 
		m_tot       = M_total1 
		m_tot2      = 0 
	
		m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total,Mlg  = ClassIest(DV2,Ve,m_tot,TWR,m_upper,m_tot2)
	
		m_tot       = M_total1
		m_upper     = 0 
		m_tot2      = M_total 
	
		m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1,Mlg1 = ClassIest(DV1,Ve,m_tot,TWR,m_upper,m_tot2)
		
		m_upper = M_total1 
		m_tot   = M_total1 
	
		m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg  = ClassIest(DV2,Ve,m_tot,TWR,m_upper,m_tot2)
			  
		while M_total_new > M_total + 0.001 : 
			M_total     = M_total_new 
			m_tot       = M_total1 
			m_tot2      = M_total_new 
			m_upper     = 0 
			
			m_frac1,m_prop1,m_RCS1,Fvac1,m_eng1,m_tank1,m_thr_str1,m_stage1,M_total1,Mlg1= ClassIest(DV1,Ve,m_tot,TWR,m_upper,m_tot2)
		
			m_upper = M_total1 
			m_tot   = M_total1 
			
			m_frac,m_prop,m_RCS,Fvac,m_eng,m_tank,m_thr_str,m_stage,M_total_new,Mlg  = ClassIest(DV2,Ve,m_tot,TWR,m_upper,m_tot2)
			   
		M_dry       = m_eng + m_thr_str + m_tank 
	
		M_dry1      = m_eng1 + m_thr_str1 + m_tank1 + Mlg1
	
		M_dry_tot   = M_dry1 + M_dry + Mass_caps + Mlg
		
		m_prop_tot  = m_prop1+m_prop
		
		return [M_dry_tot,m_prop_tot]
	
	def takeoff_wing_sizing_shuttle_like(M_takeoff,takeoff_mach,takeoff_cl):
		
		rho     = 0.02 #Density at surface on Mars
		vel     = 240 * takeoff_mach #takeoff speed as function of speed of sound (240) in m/s
		S       = M_takeoff*3.7/(0.5*rho*vel**2*takeoff_cl)
		b       = 30.5/874.5* 8 * S #Adjusted manually from spaceshuttle planform
		taper   = 0.2 #Roughly, from spaceshuttle planform again
		c_root  = S / (0.5*b)*(1/(1+taper))
		c_tip   = taper * c_root
		
		return S, b, c_root, c_tip
	
	
	def AVID_wing_mass( M_land , b , S_exp , c_root , tc ):
		
		M_land  = M_land * 2.205
		b       = b / 0.3048

		S_exp   = S_exp / (0.3048**2) * 0.8 #90% of wing exposed
	
		c_root  = c_root / 0.3048
	
		Mwing   = 1.498 * S_exp **1.176 * 0.4536 * 0.5
	
		return Mwing
	
	if k == "SPACEPLANE":
		
		m_upper     = 0
		TWR         = 1.5
		DeltaV      = (DV1+DV2)*f # Assumed less DeltaV for landing
		m_tot       = Mass_caps
		m_wing      = 1000
	
		m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total, M_dry, m_wing = Class_I_spaceplane_est(DeltaV,Ve,m_tot,TWR,m_wing)
	
		S,b,c_root,c_tip = takeoff_wing_sizing_shuttle_like(M_dry*1.1,0.8,1.5)
	
		m_wing           = AVID_wing_mass(M_dry,b,S,c_root,0.15)
	
		m_tot            = M_total
	
		m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total_new, M_dry, m_wing = Class_I_spaceplane_est(
			DeltaV, Ve, m_tot, TWR, m_wing)
	
		while M_total_new > M_total + 0.001:
			M_total = M_total_new
	
			S, b, c_root, c_tip     = takeoff_wing_sizing_shuttle_like(M_dry*1.1, 0.8, 1.5)
			m_wing                  = AVID_wing_mass(M_dry,b,S, c_root, 0.15)
	
			m_tot                   = M_total
	
			m_frac, m_prop, m_RCS, Fvac, m_eng, m_tank, m_thr_str, M_total_new, M_dry, m_wing = Class_I_spaceplane_est(
				DeltaV, Ve, m_tot, TWR, m_wing)
			
		return [M_dry,m_prop]
	
	def SpaceElevator():
		
		climbermass             = Mass_caps
		power_specific          = 1002 #W/kg
		efficiency              = 0.03*0.59*0.82
		power_per_climbermass   = 2.4*10**6/20000
		power                   = power_per_climbermass*climbermass/efficiency*(1+MGA_se_power)
		mprop_equivalent        = power/power_specific

		maxheight               = 100000000 #m
		steps                   = 51
		areostationary_height   = 17032000 #m
		numcables               = 3
		safety_factor           = 1
		A_base                  = 0.0000105*safety_factor
		height_steps            = np.linspace(areostationary_height+10,maxheight,steps)
		m_cable                 = A_base*1400*areostationary_height
		m_cables                = []
		m_totals                = []
		m_counterweights        = []
		
		for height in height_steps:
				
			taper_ratio             = exp(3389000*1400*3.71/(2*48600000000)*((3389000/height)**3-3*(3389000/height)+2))
			m_cablesegment          = taper_ratio*A_base*1400*(maxheight/(steps-1))*numcables
			m_cable                 += m_cablesegment
			m_counterweight         = 1400*A_base*48600000000*exp((3389000**2*1400*3.71)/(2*48600000000*17032000**3)*((2*17032000**3+3389000**3)/3389000-(2*17032000**3+(height)**3)/(height)))/((3389000**2*(height))/17032000**3*(1-(17032000/(height))**3)*1400*3.71)
			m_total                 = m_cable + m_counterweight + Mass_caps
			m_cables.append(m_cable)
			m_counterweights.append(m_counterweight)
			m_totals.append(m_total)
		m_total = min(m_totals)
		maxheight = height_steps[m_totals.index(min(m_totals))]
		m_counterweight = m_counterweights[m_totals.index(min(m_totals))]
		m_cable = m_cables[m_totals.index(min(m_totals))]
		return m_total, mprop_equivalent
		
	if k == "SE" : 
		
		m_total_se, mprop_equivalent = SpaceElevator()
		return [m_total_se,mprop_equivalent]

	return mass_cons

