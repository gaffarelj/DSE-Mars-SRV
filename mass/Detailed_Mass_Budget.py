import numpy as np


def Mass_conc(DV1, DV2, Isp, concept, MGA="YES", N_crew=6, N_days=4):
	conv = 0.4535

	m = dict()
	Fvac = None
	
	MGA_data = {
		"EPS": 0.117, 
		"data_hand": 0.25, 
		"TPS": 0.25, 
		"com": 0.461, 
		"GNC": 0.215, 
		"cab": 0.375, 
		"ls": 0.225, 
		"se_power": 0.5
		}
	if MGA.lower() != "yes":
		MGA_data = dict.fromkeys(MGA_data, 0)

	
	m["fuel_cells"] = (3030 * N_crew / 7) * conv
	m["bat"] = 216
	m["EPS"] = m["fuel_cells"] + m["bat"]
	
	m["life_support"] = (2444 * N_crew / 7 + 645 * N_crew + 86.4 * N_days) * conv
	
	m["com"] = (131 * N_days / 7 + 1400 * N_crew / 7) * conv
	m["data_hand"] = (302 + 828 * N_days / 7 + 1010 * N_crew / 7) * conv 
	m["GNC"] = (242 + 108 * N_days / 7 + 617 * N_crew / 7) * conv 
	
	m["av"] = m["GNC"] * (1 + MGA_data["GNC"]) + m["data_hand"] * (1 + MGA_data["data_hand"]) + m["com"] * (1 + MGA_data["com"])
	
	m["TPS"] = 2500
	
	m["cabin"] = (28.31 * (39.66 * (N_crew * N_days) ** 1.002) ** 0.6916) * conv
	
	m["pl"] = 1200
	m["payl_cont"] = 0.7 * m["pl"]
	
	g0 = 9.80665
	Ve = Isp * g0
	TWR = 1.5
	m["caps"] = m["EPS"] * (1 + MGA_data["EPS"]) + m["life_support"] * (1 + MGA_data["ls"]) + \
		m["av"] + m["cabin"] * (1 + MGA_data["cab"]) + m["pl"] + m["payl_cont"] + m["TPS"] * (1 + MGA_data["TPS"])

	def m_frac_prop(DeltaV, m_tot):
		m["frac"] = np.exp(DeltaV / Ve)    
		m["prop"] = m["frac"] * m_tot / (1 + m["frac"])

	def m_propuls(m_tot):
		Fvac = TWR * 3.7 * m_tot
		m["RCS"] = 0.0126 * m_tot
		m["eng"] = 0.00514 * Fvac ** 0.92068
		m["thr_str"] = 1.949 * 10 ** (-3) * (Fvac / 4.448) ** 1.0687 * 0.453

	def m_tank(rho_f=423, rho_ox=1140, F_o=3.8, MEOP=3e6):
		Mf = m["prop"] / (1 + 1 / F_o)
		Vf = Mf / rho_f
		Mox = m["prop"] - Mf
		Vox = Mox / rho_ox
		m["tank"] = (Vf + Vox) * MEOP / (6.43 * 10 ** 4)
	
	def ClassIest(DeltaV, m_tot, m_upper=0, m_tot2=0):
		klg = 0.033
		m_frac_prop(DeltaV, m_tot)
		m["stage"] = 0.001148 * m_upper
	
		m_tank()

		m_propuls(m_tot) if m["stage"] == 0 else m_propuls(m_tot2)

		m["dry"] = m["eng"] + m["thr_str"] + m["tank"] + m["caps"]
		m["lg"] = klg * m["dry"] * 1.15

		m["total"] = m["prop"] + m["tank"] + m["thr_str"] + m["eng"] + m["stage"] + m["RCS"] + m["lg"]
		m["total"] += m["caps"] if m["stage"] == 0 else m_tot
	
		return m["frac"], m["prop"], m["RCS"], Fvac, m["eng"], m["tank"], m["thr_str"], m["stage"], m["total"], m["lg"]
	
	def Class_I_spaceplane_est(DeltaV, m_tot, F_o=3.8, MEOP=3e6):
		m_frac_prop(DeltaV, m_tot)
		m_propuls(m_tot)
		m["landinggear"] = 0.010784 * ((m_tot - m["prop"] - m["RCS"]) * 0.453) ** 1.0861 * 0.453
		
		m_tank()

		m_propuls(m_tot)

		m["dry"] = m["caps"] + m["tank"] + m["eng"] + m["thr_str"] + m["wing"] + m["landinggear"]
		m["total"] = m["caps"] + m["prop"] + m["tank"] + m["thr_str"] + m["eng"] + m["RCS"] + m["wing"] + m["landinggear"]

		return m["frac"], m["prop"], m["RCS"], Fvac, m["eng"], m["tank"], m["thr_str"], m["total"], m["dry"], m["wing"]
	
	def takeoff_wing_sizing_shuttle_like(M_takeoff, takeoff_mach=0.8, takeoff_cl=1.5):
		rho = 0.02					# Air density at surface on Mars
		vel = 240 * takeoff_mach	# Takeoff speed as function of speed of sound (240) in m/s
		S = M_takeoff * 3.7 / (0.5 * rho * vel ** 2 * takeoff_cl)
		b = 30.5 / 874.5 * 8 * S	# Adjusted manually from spaceshuttle planform
		taper = 0.2					# Roughly from spaceshuttle planform
		c_root = S / (b/2) * (1 / (1 + taper))
		c_tip = taper * c_root
		return S, b, c_root, c_tip
	
	def AVID_wing_mass(m_land, b , S_exp , c_root , tc):
		m_land = m_land * 2.205
		b = b / 0.3048
		S_exp = S_exp / (0.3048 ** 2) * 0.8 # 80% of wing exposed
		c_root = c_root / 0.3048
		m["wing"] = 1.498 * S_exp ** 1.176 * 0.4536 * 0.5
	
	if concept.lower() == "ssto":
		DeltaV = DV1 + DV2
		
		ClassIest(DeltaV, m["caps"])
		m_old = m["total"]
		m_new = m_old + 1

		while m_new > m_old + 0.001 : 
			m_old = m_new
			ClassIest(DeltaV, m_old)
			m_new = m["total"]
			 
		m["dry"] = m["eng"] + m["thr_str"] + m["tank"] + m["caps"] + m["lg"]
		
		return m
		
	elif concept.lower() == "2_stage":
		# 1st stage
		ClassIest(DV1, m["caps"])
		m_stage1 = m["total"]
		# 2nd stage
		ClassIest(DV2, m["total"], m["total"])
		m_old = m["total"]
		m_new = m_old + 1
		while m_new > m_old + 0.001 : 
			m_old = m_new 
			m_tot = m_stage1 
			
			# 1st stage
			ClassIest(DV1, m_tot, m_tot2=m_new)
			m_stage1 = m["total"]
			m_dry1 = m["eng"] + m["thr_str"] + m["tank"] + m["lg"]
			m_prop1 = m["prop"]

			# 2nd stage
			ClassIest(DV2, m["total"], m_upper=m["total"], m_tot2=m_new)
			m_new = m["total"]
			   
		m_dry = m["eng"] + m["thr_str"] + m["tank"]
		m_dry_tot = m_dry1 + m_dry + m["caps"] + m["lg"]
		m_prop_tot = m_prop1 + m["prop"]
		
		return [m_dry_tot, m_prop_tot]

	elif concept.lower() == "spaceplane":
		DeltaV = DV1 + DV2
		m["wing"] = 1000
	
		Class_I_spaceplane_est(DeltaV, m["caps"])
		m_old = m["total"]
		S, b, c_root, c_tip = takeoff_wing_sizing_shuttle_like(m["dry"] * 1.1)
		AVID_wing_mass(m["dry"], b, S, c_root, 0.15)
		
		m_new = m_old + 1
		while m_new > m_old + 0.001:
			m_old = m_new
	
			S, b, c_root, c_tip = takeoff_wing_sizing_shuttle_like(m["dry"] * 1.1)
			AVID_wing_mass(m["dry"], b, S, c_root, 0.15)
	
			Class_I_spaceplane_est(DeltaV, m["total"])
			m_new = m["total"]
			
		return m
		
	elif concept.lower() == "se" :
		climbermass = m["caps"]
		power_specific = 1002 #W/kg
		efficiency = 0.03 * 0.59 * 0.82
		power_per_climbermass = 2.4 * 10 ** 6 / 20000
		power = power_per_climbermass * climbermass / efficiency * (1 + MGA_data["se_power"])
		mprop_equivalent = power / power_specific

		maxheight = 1e8 #m
		steps = 51
		areostationary_height = 17.032e6 #m
		numcables = 3
		safety_factor = 1
		A_base = 0.0000105 * safety_factor
		height_steps = np.linspace(areostationary_height + 10, maxheight, steps)
		m_cable = A_base * 1400 * areostationary_height
		m_cables = []
		m_totals = []
		m_counterweights = []
		modulus = 48.6e9
		
		for height in height_steps:
			taper_ratio = np.exp(3389e3 * 1400 * 3.71 / (2 * modulus) * ((3389e3 / height) ** 3 - 3 * (3389e3 / height) + 2))
			m_cablesegment = taper_ratio * A_base * 1400 * (maxheight / (steps - 1)) * numcables
			m_cable += m_cablesegment
			m_counterweight = 1400 * A_base * modulus * np.exp((3389e3 ** 2 * 1400 * 3.71) / (2 * modulus * areostationary_height ** 3) 
				* ((2 * areostationary_height ** 3 + 3389000 ** 3) / 3389000 - (2 * areostationary_height ** 3 + (height) ** 3) / (height))) / ((3389000 ** 2 
				* (height)) / areostationary_height ** 3 * (1 - (areostationary_height / (height)) ** 3) * 1400 * 3.71)
			m_total = m_cable + m_counterweight + m["caps"]
			m_cables.append(m_cable)
			m_counterweights.append(m_counterweight)
			m_totals.append(m_total)
		m_total = min(m_totals)
		maxheight = height_steps[m_totals.index(min(m_totals))]
		m_counterweight = m_counterweights[m_totals.index(min(m_totals))]
		m_cable = m_cables[m_totals.index(min(m_totals))]
		return m_total, mprop_equivalent
	
	raise Exception("Undefined concept name")