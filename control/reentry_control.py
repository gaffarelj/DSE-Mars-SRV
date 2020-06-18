import sys
sys.path.append('../astrodynamics')
# import mars_standard_atmosphere as MSA
# from reentry_footprint import flight, time, dt, mean_radius, mars
import disturbances as dist
import actuator_properties as act
import pitching_moment as pm

import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib import pyplot as plt


def reentry_control(I,cg,tburn,thrust_level,Isp):
    margin = 2.

    plotting = False
    printing = True
    #=====================================================================================================================================================================================================================
    #Vehicle onstants
    #=====================================================================================================================================================================================================================
    Ix = I[0]
    Iy = I[1]
    Iz = I[2]

    g   = act.g
    thrust_level = 451.85155228408263

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
    R_mars   = 3389.5*10**3

    #=====================================================================================================================================================================================================================
    #Disturbance
    #=====================================================================================================================================================================================================================
    angle = 2. * np.pi /180
    omega = velocity / R_mars
    S, cp = dist.Drag_surface_cp(angle)
    drag = dist.Drag_force(q,cd,S)
    Td = dist.aerodynamic_disturbance(cp,cg,drag,angle)
    Tgx = dist.gravitygradient_disturbance(Iy,Iz,omega,angle)
    Tgy = dist.gravitygradient_disturbance(Ix,Iz,omega,angle)
    Tsp= dist.solarpressure_disturbance(angle,cg)
    Tm = dist.magnetic_disturbance(R_mars)
    T_dist = Tgx + Tsp + Tm


    #=====================================================================================================================================================================================================================
    #Function to calculate total impulse required for slew maneuver
    #=====================================================================================================================================================================================================================
    def slew_landing(thrust,alpha0,S,cp,cd,q,T_dist,cg,I,t0,t_end,Isp,g,margin):
        slew_angle_tot = 180 * np.pi / 180 - alpha0
        thrust = thrust
        slew_duration = t_end - t0
        spin_rate_avg = slew_angle_tot / slew_duration
        t             = 0
        dt            = 0.1
        slew_angle    = 0
        spin_rate     = 0
        impulse       = 0
        mp            = 0

        while slew_angle < slew_angle_tot/2 and t < (t_end-t0)/2:
            alpha         = alpha0 + slew_angle
            # if sum(spin_rates)/len(spin_rates) > spin_rate_avg:
            #     thrust = RCS_thrust * 0
            # else:
            #     thrust = RCS_thrust
            RCS_torque    = act.RCS_thrust_to_torque(thrust,"z",cg)

            net_torque    = RCS_torque - T_dist
            # print(T_dist,RCS_torque)
            spin_acc      = (net_torque) / I
            spin_rate    += spin_acc * dt

            slew_angle   += spin_rate * dt
            t            += dt
            impulse      += thrust * dt
            mp           += (thrust * dt) / (Isp*g)

        #Rotation successfully completed or not
        if slew_angle >= slew_angle_tot/2:
            success = 'yes'
        else:
            success = 'no'
        mp = mp * 2
        return impulse,t, mp, success,RCS_torque

    #=====================================================================================================================================================================================================================
    #All possible rotations with corresponding impulse and time
    #=====================================================================================================================================================================================================================

    rotation_values = []
    for t0 in range(int(t_end-100),int(t_end-10)):

        impulse, slew_time, mp, success, RCS_torque = slew_landing(thrust_level,alpha,S,cp,cd,q,T_dist,cg,Iy,t0,t_end,Isp,g,margin)
        mp     = mp * margin
        if success == 'yes':
            # print([thrust,impulse,slew_time,mp])
            rotation_values.append([impulse,slew_time,mp,RCS_torque])


    rotation_values = np.array(rotation_values)
    mp     = rotation_values[-1,2]
    thrust = thrust_level
    t      = rotation_values[-1,1]
    torque = rotation_values[-1,-1]

    #=====================================================================================================================================================================================================================
    #Errors
    #=====================================================================================================================================================================================================================
    #========================================================================
    # Thruster misalignment
    #========================================================================
    error_angle = 2 * np.pi / 180
    n_torquethrusters = 18.
    T_error_y, T_error_x, T_error_z = act.thrust_error(thrust,cg,error_angle)

    RCS_error_x  = margin * act.RCS_torque_to_thrust(T_error_x,'y',cg,'error_bottom')
    RCS_error_y  = margin * act.RCS_torque_to_thrust(T_error_y,'x',cg,'error_bottom')
    RCS_error_z  = margin * act.RCS_torque_to_thrust(T_error_z,'z',cg,'error_bottom')

    RCS_error    = max([RCS_error_x,RCS_error_y,RCS_error_z])


    #=====================================================================================================================================================================================================================
    #Aerodynamic pitch control
    #=====================================================================================================================================================================================================================
    t_pitch = pm.time

    pitch_moments = pm.pitching_moment
    pitch_mp_tot = []
    pitch_mp     = 0
    pitch_thrust_max = 4 *act.RCS_torque_to_thrust(abs(min(pitch_moments)),'z',cg,'normal')
    for pitch_moment in pitch_moments:
        pitch_thrust = 4 * act.RCS_torque_to_thrust(abs(pitch_moment),'z',cg,'normal')
        pitch_mp = act.RCSpropellant(pitch_thrust,0.1,Isp)
        pitch_mp_tot.append(pitch_mp)

    #=====================================================================================================================================================================================================================
    #Total
    #=====================================================================================================================================================================================================================
    T_error      = [T_error_x, T_error_y, T_error_z]
    RCS_error    = max([RCS_error_x,RCS_error_y,RCS_error_z])
    mp_error     = n_torquethrusters * act.RCSpropellant(RCS_error,t,Isp)

    T_rot        = [torque, 0., 0.]
    RCS_rot      = thrust
    mp_rot       = np.sum(mp)
    t_rot        = t

    T_pitch      = [abs(min(pm.pitching_moment)), 0., 0.]
    RCS_pitch    = pitch_thrust_max
    mp_pitch     = sum(np.array(pitch_mp_tot))
    t_pitch      = t_pitch[-1]

    mp_total     = mp + mp_error + pitch_mp

    if printing == True:
        print('==========================================================')
        print('==========================================================')
        print('REENTRY')
        print('')
        print('ROTATION')
        print('Total torque (x,y,z)           : ', T_rot)
        print('Thrust per engine (x,y,z)      : ', RCS_rot)
        print('propellant needed              : ', mp_rot)
        print('time                           : ', t_rot)
        print('REDUNDANCY')
        print('Disturbance torque             : ', T_error)
        print('redundancy thrust per engine   : ', RCS_error)
        print('redundacy propellant           : ', mp_error)
        print('')
        print('PITCH CONTROL')
        print('Max pitching moment            : ', T_pitch)
        print('Max thrust                     : ', RCS_pitch)
        print('Propellant needed              : ', mp_pitch)
        print('time                           : ', t_pitch)
        print('')
        print('Total propellant mass          :', mp_total)
        print('==========================================================')
        print('==========================================================')

    #=====================================================================================================================================================================================================================
    #Plot
    #=====================================================================================================================================================================================================================
    if plotting == True:
        plt.figure()
        plt.plot(t_pitch,pitch_mp_tot,color="navy")
        plt.grid(color="gainsboro")
        plt.title("Time vs Propellant mass")
        plt.xlabel("Time [s]")
        plt.ylabel("Propellant mass [kg]")
        plt.show()

    return mp_total, mp_rot, mp_pitch, RCS_rot+RCS_error, RCS_pitch, t_rot
