import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act

import numpy as np
from matplotlib import pyplot as plt

def ascent_control(tm,am,I,cg_orbit,tburn,RCS_thrust_rot,Isp):
    margin = 2.

    printing = True
    plotting = False
    #=====================================================================================================================================================================================================================
    #Node properties
    #=====================================================================================================================================================================================================================
    mu     = 0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
    R      = 3389.5*10**3                  #[m] volumetric mean radius
    h_node = 500*10**3
    h_phasing=609.74*10**3
    g_earth= 9.81             #[m]

    period =  2*np.pi* np.sqrt((R+h_node+h_phasing) ** 3 / mu)
    omega  = (2 * np.pi)/period
    V      = omega * (R+h_node+h_phasing)

    #=====================================================================================================================================================================================================================
    #Initial slew and Delta V values
    #=====================================================================================================================================================================================================================
    #Assume spin acceleration/deceleration of 5%, coast time of 90%
    angle = (180. - 26.7) * np.pi / 180
    slew_limit = 0.5 * period #s

    #=====================================================================================================================================================================================================================
    #Vehicle constants
    #=====================================================================================================================================================================================================================
    length = act.length_body
    width  = act.body_radius * 2
    Ix = I[0]
    Iy = I[1]
    Iz = I[2]
    m  = 53056.28092

    #=====================================================================================================================================================================================================================
    #Function to compute thrust from slew maneuver
    #=====================================================================================================================================================================================================================
    def slew_ascent(tburn,slew_angle,I,cg,Isp,RCS_thrust):
        torque = act.RCS_thrust_to_torque(RCS_thrust,'z',cg)
        spin_acc   = torque / I
        spin_rate  = spin_acc * tburn
        slew_time = slew_angle / spin_rate
        RCS_impulse = 4 * RCS_thrust * tburn * 2
        Mp = RCS_impulse / (Isp * 9.80665)
        return torque, RCS_thrust, RCS_impulse, Mp,slew_time


    #========================================================================
    # Slew maneuver
    #========================================================================
    RCS_torque, RCS_thrust, RCS_impulse, Mp,slew_time = slew_ascent(tburn,angle,Iy,cg_orbit,Isp,RCS_thrust_rot)
    #========================================================================
    # Aerodynamic pitch control
    #========================================================================
    dt = 0.1
    pitch_mp = 0
    pitch_thrust_max = act.RCS_torque_to_thrust(abs(min(am)),'z',cg_orbit,'normal')

    pitch_mp_tot = []
    pitch_mp = 0
    for aero_moment in am:
        pitch_thrust = 2*act.RCS_torque_to_thrust(aero_moment,'z',cg_orbit,'normal')
        pitch_mp    += 4 * act.RCSpropellant(abs(pitch_thrust),dt,Isp)
        pitch_mp_tot.append(pitch_mp)
    pitch_mp_tot = np.array(pitch_mp_tot)


    #=====================================================================================================================================================================================================================
    #Errors
    #=====================================================================================================================================================================================================================
    #========================================================================
    # Thruster misalignment
    #========================================================================
    error_angle = 2 * np.pi / 180

    T_error_y, T_error_x, T_error_z = act.thrust_error(RCS_thrust,cg_orbit,error_angle)
    T_error_pitch_y, T_error_pitch_x, T_error_pitch_z = act.thrust_error(pitch_thrust_max,cg_orbit,error_angle)
    #========================================================================
    # Disturbances
    #========================================================================
    theta = 2. * np.pi / 180

    Tgx  = dist.gravitygradient_disturbance(Iy,Iz,omega,theta)
    Tgy  = dist.gravitygradient_disturbance(Ix,Iz,omega,theta)
    Tsp = dist.solarpressure_disturbance(theta,cg_orbit)
    Tm  = dist.magnetic_disturbance(R)

    #========================================================================
    # Total
    #========================================================================
    T_error_x += Tgx + Tsp + Tm
    T_error_y += Tgy
    RCS_error_x  = margin * act.RCS_torque_to_thrust(T_error_x,'y',cg_orbit,'error_bottom')
    RCS_error_y  = margin * act.RCS_torque_to_thrust(T_error_y,'x',cg_orbit,'error_bottom')
    RCS_error_z  = margin * act.RCS_torque_to_thrust(T_error_z,'z',cg_orbit,'error_bottom')

    RCS_error    = max([RCS_error_x,RCS_error_y,RCS_error_z])
    mp_error     = act.RCSpropellant(RCS_error,tburn,Isp)

    T_error_pitch_x += Tgx + Tsp + Tm
    T_error_pitch_y += Tgy
    RCS_error_pitch_x  = margin * act.RCS_torque_to_thrust(T_error_pitch_x,'y',cg_orbit,'error_bottom')
    RCS_error_pitch_y  = margin * act.RCS_torque_to_thrust(T_error_pitch_y,'x',cg_orbit,'error_bottom')
    RCS_error_pitch_z  = margin * act.RCS_torque_to_thrust(T_error_pitch_z,'z',cg_orbit,'error_bottom')

    RCS_pitch_error  = max([RCS_error_pitch_x,RCS_error_pitch_y,RCS_error_pitch_z])
    mp_pitch_error     = 18 * act.RCSpropellant(RCS_pitch_error,tm[-1],Isp)
    #=====================================================================================================================================================================================================================
    #Total
    #=====================================================================================================================================================================================================================
    T_rot   = [RCS_torque, 0., 0.]
    RCS_rot = RCS_thrust
    mp_rot  = Mp
    t_rot   = slew_time

    T_error_rot     = [T_error_x, T_error_y, T_error_z]
    RCS_error_rot   = max([RCS_error_x,RCS_error_y,RCS_error_z])
    mp_error_rot    = mp_error

    T_error_pitch   = [T_error_pitch_x, T_error_pitch_y, T_error_pitch_z]
    RCS_error_pitch = max([RCS_error_pitch_x,RCS_error_pitch_y,RCS_error_pitch_z])
    mp_error_pitch  = mp_pitch_error

    mp_error_tot    = mp_error_rot + mp_error_pitch

    T_pitch         = [abs(min(am)), 0., 0.]
    RCS_pitch       = [pitch_thrust_max, 0., 0.]
    mp_pitch        = pitch_mp_tot[-1]
    t_pitch         = tm[-1]

    mp_tot          = mp_rot + mp_error_tot + mp_pitch
    if printing == True:
        print('==========================================================')
        print('==========================================================')
        print('ASCENT')
        print('')
        print('ROTATION')
        print('Torque                                      : ', T_rot)
        print('Thrust per engine                           : ', RCS_rot)
        print('propellant needed                           : ', mp_rot)
        print('time                                        : ', t_rot)
        print('time limit                                  : ', slew_limit)
        print('')
        print('PITCH CONTROL')
        print('Max pitch moment                            : ', T_pitch)
        print('Max thrust                                  : ', RCS_pitch)
        print('propellant needed                           : ', mp_pitch)
        print('time                                        : ', tm[-1])
        print('')
        print('REDUNDANCY')
        print('Disturbance torque (rotation)               : ', T_error_rot)
        print('Redundancy thrust per engine (rotation)     : ', RCS_error_rot)
        print('Redundancy propellant (rotation)            : ', mp_error_rot)
        print('Disturbance torque (pitch control)          : ', T_error_pitch)
        print('Redundancy thrust per engine (pitch control): ', RCS_error_pitch)
        print('Redundancy propellant (pitch control)       : ', mp_error_pitch)
        print('')
        print('Total propellant needed                     : ',mp_tot)
        print('==========================================================')
        print('==========================================================')
    #=====================================================================================================================================================================================================================
    #Plot
    #=====================================================================================================================================================================================================================
    if plotting == True:
        plt.figure()
        plt.plot(tm[1:],pitch_mp_tot,color="navy")
        plt.grid(color="gainsboro")
        plt.title("Time vs Propellant mass")
        plt.xlabel("Time [s]")
        plt.ylabel("Propellant mass [kg]")
        plt.show()

    return mp_tot, mp_rot+mp_error_rot, mp_pitch+mp_error_pitch, RCS_rot+RCS_error_rot, RCS_pitch+RCS_error_pitch,t_rot
