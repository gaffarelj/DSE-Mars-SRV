import numpy as np
from matplotlib import pyplot as plt
import sys
import LAS_ascent as LASA
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Aero'))
import aero_calcs as ac
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'astrodynamics'))
import astro_tools
import astro_tools_nonfull as AT

mean_radius = 3389500							#[m]
reentry_altitude = 80000                        #[m]
v_insertion_orbit = 3272.46591                  #[m/s] 
insertion_orbit_a = mean_radius + 500000        #[m]
insertion_orbit_p = mean_radius - 300000        #[m]

vehicle_mass = 38000                            #[kg]
S = 350                                         #[m^2]
MOI =  [391281.1045, 22714482.967, 22714482.967]    # Ixx, Iyy, Izz
dt = 0.1 

mars = astro_tools.Planet()

state = np.zeros(12)
state[0] = mars.reentry_velocity(reentry_altitude, insertion_orbit_a, insertion_orbit_p)                        # velocity
state[1] = mars.reentry_angle(reentry_altitude, insertion_orbit_p, insertion_orbit_a)                           # flight path angle (1.5 arcsecs (2 sigma))
state[2] = np.radians(42.5)                                                                                     # heading angle (1.5 arcsecs (2 sigma))
state[3] = reentry_altitude + mean_radius                                                                       # radius 
state[4] = np.radians(-1.168)                                                                                   # lognitude 25.5 -1.168 - 0.141
state[5] = np.radians(24.322)                                                                                   # latitude 42.5 24.322 - 0.077
state[6] = 0																									# rollrate 
state[7] = 0																									# pitchrate
state[8] = 0																									# yawrate
state[9] = -np.radians(60) 																						# angle of attack defined positive downwards 
state[10] = 0																									# sideslip angle  
state[11] = np.radians(0)																									# bank angle 	

run = False
if run:
	motion = astro_tools.Motion(state, MOI, S, vehicle_mass, ac.H_aerodynamics_coefficients, mars) #, True, 230e3, 733.45)
	flight, time = motion.forward_euler(dt)


	#astro_tools.plot_single(time, motion.heatflux, 'Time [s]', 'Heat Flux [W/m^2]')
	heat_ok_i = np.where(np.array(motion.heatflux) <= 2000)[0][0]
	# time, velocity, altitude, gamma, AoA
	t, V0, h0, gamma, AoA = time[heat_ok_i], flight[heat_ok_i][0], flight[heat_ok_i][3]-mean_radius, flight[heat_ok_i][1], flight[heat_ok_i][9]
	print(t, V, h, gamma, AoA)
else:
	t, V0, h0, gamma, AoA = 601.5, 1009.46, 22898, 0, -0.5

m_fuel_tot = 2336
m_fuel_abort = 1840
Isp_abort_opt = 221.00
Isp_abort_vac = 259.16
m = 14200
vehicle_mass = m - m_fuel_abort
dV_abort = LASA.dV(Isp_abort_opt, m, m_fuel_abort)
V0 -= dV_abort
q_s = [460, 500]
chute_ts = [9999, 9999]
for i, q_chute in enumerate(q_s):
	print(f"Finding deploy time for parachute {i+1} (at q={q_chute} Pa).")
	chutes = LASA.def_chutes(chute_ts)
	t, a, v, q, h, angle, m = LASA.run_motion(V0, gamma*180/np.pi, h0, chutes, vehicle_mass=vehicle_mass, aoa=AoA*180/np.pi)
	Vy = v * np.sin(angle)
	# Index for which Vy = 0
	v_0_i = 0#np.where(Vy <= 0)[0][0]
	# Time at which q = q_chute
	chute_ts[i] = t[np.where(q[v_0_i:] <= q_chute)[0][0] + v_0_i]
	
chute_ts[1] = chute_ts[0] + 24
chutes = LASA.def_chutes(chute_ts)
print("Parachutes deployment times are at ", ", ".join([str(round(t, 2)) for t in chute_ts]), "[s]")
t, a, v, q, h, angle, final_m = LASA.run_motion(V0, gamma*180/np.pi, h0, aoa=AoA*180/np.pi, chutes=chutes, vehicle_mass=vehicle_mass, print_deploy=True)

dV_land = LASA.dV(Isp_abort_vac, final_m, m_fuel_tot - m_fuel_abort)

print(f"Maximum acceleration was {round(min(a)/9.81, 2)} g0.")
print(f"Landing with velocity of {round(v[-1], 2)} m/s, with {round(dV_land, 2)} m/s for propulsive landing.")

AT.plot_dual(t, q, v, 'Time [s]', 'Dyn.  Pressure [Pa]', 'Velocity [m/s]')
AT.plot_dual(t, v, [ht / 1000 for ht in h], 'Time [s]', 'Velocity [m/s]', 'Altitude [km]')
AT.plot_dual(t, a, q, 'Time [s]', 'Acceleration [m/s$^2$]', 'Dynamic Pressure [Pa]')