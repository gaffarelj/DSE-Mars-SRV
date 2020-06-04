import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT
from astrodynamics import BAT


def al(n):
	try:
		return n[0]
	except IndexError:
		return n

def run_ascent(max_M=7, min_q=100):
	ascent = BAT.ascent_sim()[4]
	# Get index when M = 8
	Mach_T_i = np.where(np.array(ascent["Mach"]) >= max_M)[0][0]
	# Get index when q = 50 Pa (after max reached)
	q_max_i = np.where(ascent["q"] == max(ascent["q"]))[0][0]
	q_T_i = np.where(ascent["q"][q_max_i:] <= min_q)[0][0] + q_max_i
	# Take minimum
	LAS_T_i = min(q_T_i, Mach_T_i)
	# Get corresponding time, and data
	LAS_T = round(ascent["time"][LAS_T_i], 2)
	Mach_T = round(al(ascent["Mach"][LAS_T_i]), 2)
	q_T = round(ascent["q"][LAS_T_i], 2)
	V_T = round(ascent["V"][LAS_T_i], 2)
	h_T = round(al(ascent["altitude"][LAS_T_i]), 2)
	# Print results
	gamma = ascent["gamma"][-1]
	print(f"Last abort to surface at t={LAS_T} s, at M={Mach_T} and q={q_T} Pa")
	print(f"\t- at a velocity of {V_T} m/s and an altitude of {h_T} m")
	return V_T, gamma, h_T

def run_motion(V0, gamma, h0, chutes=[], print_deploy=False, prop_reentry=[]):
	mars = AT.Planet()
	state = np.zeros(6)
	state[0] = V0						# velocity
	state[1] = np.radians(gamma)		# flight path angle
	state[2] = np.radians(gamma)		# heading angle
	state[3] = mars.r + h0				# radius 
	state[4] = np.radians(-27.088)      # longitude 25.5
	state[5] = np.radians(4.51)         # latitude 42.5

	roll_angle = np.radians(0)
	aoa = np.radians(-5)

	cd = 1.2
	cl = 0.23
	S = 11.9

	motion = AT.Motion(state, roll_angle, aoa, S, vehicle_mass, cl, cd, mars, chutes, print_deploy, prop_reentry)
	flight, time = motion.forward_euler(0.1)
	return time, np.array(motion.a_s), flight[:,0], np.array(motion.q_s), flight[:,3] - mars.r, flight[:,1], motion.mass

def def_chutes(times):
	#ballute = AT.pc(0.35, 12, 50, deploy_time=times[0], n=1, name="ballute")
	drogue = AT.pc(0.35, 10, 20, deploy_time=times[0], n=5, name="drogue")
	main = AT.pc(0.55, 25, 30, deploy_time=times[1], n=3, name="main")
	return [drogue, main]

def dV(Isp, m0, mf):
	return 9.81*Isp*np.log(m0/(m0-mf))


#V_T, gamma, h_T = run_ascent()
#print(V_T, gamma, h_T)
V_T, gamma, h_T = 1427.87, 48.02, 43497.6

m_fuel_tot = 2000
m_fuel_abort = 1176.84
Isp_abort_opt = 221.2
Isp_abort_vac = 259.38
m = 12500
vehicle_mass = m - m_fuel_abort
dV_abort = dV(Isp_abort_opt, m, m_fuel_abort)

# Dynamic pressures at which the parachutes will be deployed
q_s = [20, 350]
chute_ts = [9999, 9999]
for i, q_chute in enumerate(q_s):
	print(f"Finding deploy time for parachute {i+1} (at q={q_chute} Pa).")
	chutes = def_chutes(chute_ts)
	t, a, v, q, h, angle, m = run_motion(V_T + dV_abort, gamma, h_T, chutes)
	Vy = v*np.sin(angle)
	# Index for which Vy = 0
	v_0_i = np.where(Vy <= 0)[0][0]
	# Time at which q = q_chute
	chute_ts[i] = t[np.where(q[v_0_i:] >= q_chute)[0][0] + v_0_i]
	
chutes = def_chutes(chute_ts)
print("Parachutes deployment times are at ", ", ".join([str(round(t, 2)) for t in chute_ts]), "[s]")
t, a, v, q, h, angle, final_m = run_motion(V_T + dV_abort, gamma, h_T, chutes, print_deploy=True)

dV_land = dV(Isp_abort_vac, final_m, m_fuel_tot-m_fuel_abort)

print(f"Maximum acceleration was {round(min(a)/9.81, 2)} g0.")
print(f"Landing with velocity of {round(v[-1], 2)} m/s, with {round(dV_land, 2)} m/s for propulsive landing.")

#AT.plot_dual(t, q, v, 'Time [s]', 'Dyn. Pressure [Pa]', 'Velocity [m/s]')
AT.plot_dual(t, q, [ht/1000 for ht in h], 'Time [s]', 'Dyn. Pressure [Pa]', 'Altitude [km]')
AT.plot_dual(t, [at/9.81 for at in a], v, 'Time [s]', 'Acceleration [g0]', 'Velocity [m/s]')