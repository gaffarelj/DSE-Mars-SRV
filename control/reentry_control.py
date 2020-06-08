import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
# from reentry_footprint import flight, time, dt, mean_radius, mars
import disturbances as dist
import actuator_properties as act

import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib import pyplot as plt

#=====================================================================================================================================================================================================================
#Vehicle onstants
#=====================================================================================================================================================================================================================
Iy = act.Iy
Iz = act.Iz
cg = act.z_cg_empty
# length = act.length
# width = act. width
Isp = act.Isp
g   = act.g
thrust_levels = np.arange(1,50,1)

#=====================================================================================================================================================================================================================
#Flight profile
#=====================================================================================================================================================================================================================
# height = flight[:,3]-mean_radius
# velocity = flight[:,0]

cd = 2.5
t_end = 738.5 #s
q  = 800.  #900. #Pa
alpha = 45 *np.pi / 180
velocity = 500. #m/s

#=====================================================================================================================================================================================================================
#Drag disturbance
#=====================================================================================================================================================================================================================

S, cp = dist.Drag_surface_cp(alpha)
drag = dist.Drag_force(q,cd,S)
Td = dist.aerodynamic_disturbance(cp,cg,drag,alpha)



#=====================================================================================================================================================================================================================
#Function to calculate total impulse required for slew maneuver
#=====================================================================================================================================================================================================================
def slew_landing(thrust,alpha0,S,cp,cd,q,Td0,cg,I,t0,t_end):
    slew_angle_tot = 180 * np.pi / 180 - alpha0

    slew_duration = t_end - t0
    spin_rate_avg = slew_angle_tot / slew_duration
    t             = 0
    dt            = 0.1
    slew_angle    = 0
    spin_rate     = 0
    impulse       = 0
    spin_rates    = [0]

    while slew_angle < slew_angle_tot/2 and t < (t_end-t0)/2:
        alpha         = alpha0 + slew_angle
        # if sum(spin_rates)/len(spin_rates) > spin_rate_avg:
        #     thrust = RCS_thrust * 0
        # else:
        #     thrust = RCS_thrust
        RCS_torque    = act.RCS_thrust_to_torque(thrust,"y",cg)

        S_new, cp_new = dist.Drag_surface_cp(alpha)
        drag_new      = dist.Drag_force(q,cd,S_new)
        Td_new        = dist.aerodynamic_disturbance(cp_new,cg,drag_new,alpha)
        dTd           = Td_new - Td0
        net_torque    = RCS_torque
        spin_acc      = (net_torque) / I
        spin_rate    += spin_acc * dt

        slew_angle   += spin_rate * dt
        t            += dt
        impulse      += thrust * dt
        spin_rates.append(spin_rate)
    #Rotation successfully completed or not

    if slew_angle >= slew_angle_tot/2:
        success = 'yes'
    else:
        success = 'no'

    return impulse,t,success

#=====================================================================================================================================================================================================================
#All possible rotations with corresponding impulse and time
#=====================================================================================================================================================================================================================

rotation_values = []
for t0 in range(int(t_end-200),int(t_end-50)):

    for thrust in thrust_levels:

            impulse, slew_time, success = slew_landing(thrust,alpha,S,cp,cd,q,Td,cg,Iy,t0,t_end)
            mp                          = 6 * impulse / (Isp * g)
            print(thrust,mp,success,slew_time)
            if success == 'yes':
                if thrust < 11250.:
                    rotation_values.append([thrust,impulse,slew_time,mp])


rotation_values = np.array(rotation_values)

#=====================================================================================================================================================================================================================
#Redundancy roll
#=====================================================================================================================================================================================================================
angle = 90 * np.pi / 180.
time_roll = 5.

RCS_roll_torque = act.slew(angle,time_roll,Iz)
RCS_roll_thrust = act.RCS_torque_to_thrust(RCS_roll_torque,'z',cg,'normal')
mp_roll         = act.RCSpropellant(RCS_roll_thrust*4,time_roll,Isp)
print(mp_roll,RCS_roll_thrust)
