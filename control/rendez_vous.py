import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act

import numpy as np
import sympy as sp
from matplotlib import pyplot as plt


def RV_Docking_control(I,cg,tburn_deltaV,tburn_rot,Isp,m0,error):
    margin = 2.

    #=====================================================================================================================================================================================================================
    # Switches
    #=====================================================================================================================================================================================================================
    plotting = True
    printing = True
    #=====================================================================================================================================================================================================================
    # Function to compute Thrust
    #=====================================================================================================================================================================================================================
    def vac_thrust(DeltaV,Isp,Mbegin,tb,De=0,pe=0):
        """computes vacuum thrust from DeltaV.
            Mbegin=mass at the start of the maneuver, tb=burn time, De=exit diameter of engine/thruster, pe=exhaust exit pressure of engine/thruster
        """
        Ae=np.pi/4*De*De
        thrust=(np.exp(DeltaV/(Isp*9.80665))*Mbegin-Mbegin)/(tb*np.exp(DeltaV/(Isp*9.80665)))*Isp*9.80665+Ae*(pe)
        Mprop=(np.exp(DeltaV/(Isp*9.80665))*Mbegin-Mbegin)/(np.exp(DeltaV/(Isp*9.80665)))
        return thrust, Mprop

    #=====================================================================================================================================================================================================================
    #Node properties
    #=====================================================================================================================================================================================================================

    mu     = 0.042828*10**6*(10**3)**3    #[m^3/s^2] gravitational parameter
    R      = 3389.5*10**3                  #[m] volumetric mean radius
    h_node = 500*10**3
    h_phasing=609.74*10**3
    g_earth= 9.81             #[m]

    period =  2*np.pi* np.sqrt((R+h_node) ** 3 / mu)
    period_phasing = 2*np.pi* np.sqrt((R+h_node+h_phasing) ** 3 / mu)
    omega  = (2 * np.pi)/period
    time_hohmann = np.pi*np.sqrt(((h_node+h_phasing+2*R)/2)**3*1/mu)
    V      = np.sqrt((R+h_node) / mu)

    #=====================================================================================================================================================================================================================
    #Vehicle properties
    #=====================================================================================================================================================================================================================

    #Implement Mass of Charon at t0!
    OF_ratio = 3.8
    Ix = I[0] #kg/m^2
    Iy = I[1] #kg/m^2
    Iz = I[2] #kg/m^2
    cg_orbit = cg

    #=====================================================================================================================================================================================================================
    #Vbar approach properties
    #=====================================================================================================================================================================================================================

    #Proximity operations A:
    x_A = [-1000, -250]
    t_ref_A  =  period    #s

    x0_A = x_A[0]         #m
    x1_A = x_A[1]         #m
    t_A  = t_ref_A        #s
    Vx_A   = (x1_A - x0_A) / t_A #m/s
    deltaV_A_0 = deltaV_A_1 = Vx_A * margin

    #Proximity operations B:
    x_B = [x1_A, -30]    #m
    t_ref_B  = 60          #min

    x0_B = x_B[0]          #m
    x1_B = x_B[1]          #m
    t_B  = t_ref_B * 60    #s
    Vx_B  = (x1_B - x0_B) / t_B
    deltaV_B_0 = deltaV_B_1 = Vx_B * margin
    #Docking:
    x_d  = [x1_B, -3]
    t_d_ref  =  5                              #min

    x0_d = x_d[0]                              #m
    x1_d = x_d[1]                              #m
    t_d  = (t_d_ref*((x1_d - x0_d)/10)) * 60   #s
    Vx_d  = (x1_d - x0_d) / t_d
    deltaV_d_0 = deltaV_d_1 = Vx_d * margin

    y0 = z0 = 0

    deltaV_tot = deltaV_A_0 + deltaV_A_1 + deltaV_B_0 + deltaV_B_1 + deltaV_d_0 + deltaV_d_1

    #=====================================================================================================================================================================================================================
    #Navigation measurement errors
    #=====================================================================================================================================================================================================================
    error_xm = error[0]
    error_ym = error[1]
    error_zm = error[2]

    error_xdot = error[3]
    error_ydot = error[4]
    error_zdot = error[5]

    def error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t):
        delta_x = error_xm + 6 * error_zm * (omega * t - np.sin(omega * t)) +  error_xdot * (4 / omega * np.sin(omega * t) - 3 * t) + 2 / omega * error_zdot * (1 - np.cos(omega * t))
        delta_xdot = 6 * error_zm * (omega - omega * np.cos(omega * t)) +  error_xdot * (4 * np.cos(omega * t) - 3) + 2 / omega * error_zdot * (omega *  np.sin(omega * t))
        delta_xdotdot = 6 * error_zm * (omega ** 2 * np.sin(omega * t)) +  error_xdot * (-4 * omega * np.sin(omega * t)) + 2 / omega * error_zdot * (omega ** 2 * np.cos(omega * t))

        return delta_x, delta_xdot, delta_xdotdot

    def error_y(error_ym,omega,t):
        delta_y = error_ym * np.cos(omega * t) + 1 / omega * error_ydot * np.sin(omega * t)
        delta_ydot = -omega * error_ym * np.sin(omega * t) + error_ydot * np.cos(omega * t)
        delta_ydotdot = -omega ** 2 * error_ym * np.cos(omega * t) - omega * error_ydot * np.sin(omega * t)

        return delta_y, delta_ydot, delta_ydotdot

    def error_z(error_xm,error_zm,error_zdot,omega,t):
        delta_z = error_zm * (4 - 3 * np.cos(omega * t)) + 2 / omega * error_xdot * (np.cos(omega * t) -1) + 1 / omega * error_zdot * np.sin(omega * t)
        delta_zdot=3*error_zm*omega*np.sin(omega*t)-2*error_xdot*np.sin(omega*t)+error_zdot*np.cos(omega*t)
        delta_zdotdot=3*error_zm*omega*omega*np.cos(omega*t)-2*error_xdot*omega*np.cos(omega*t)-error_zdot*omega*np.sin(omega*t)

        return delta_z, delta_zdot, delta_zdotdot


    #=====================================================================================================================================================================================================================
    # Hill equations of motion
    #=====================================================================================================================================================================================================================

    def thrust_x(x,zdot,xdotdot,omega,t):
        gamma_x = (xdotdot - 2 * omega * zdot)
        return gamma_x

    def thrust_y(y,ydotdot,omega,t):
        gamma_y = (ydotdot + omega ** 2 * y)
        return gamma_y

    def thrust_z(z,zdotdot,xdot,omega,t):
        gamma_z = (zdotdot + 2 * omega * xdot - 3 * omega ** 2 * z)
        return gamma_z

    #=====================================================================================================================================================================================================================
    # Simulation
    #=====================================================================================================================================================================================================================
    dt = 1.
    t  = 0.1
    ta = tb = td = 0.1
    tburn = tburn_deltaV
    mp = 0.
    m = m0
    m_deltaV = [m]

    delta_x, delta_xdot, delta_xdotdot = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
    delta_y, delta_ydot, delta_ydotdot = error_y(error_ym,omega,t)
    delta_z, delta_zdot, delta_zdotdot = error_z(error_xm,error_zm,error_zdot,omega,t)
    x   = x0_A + delta_x
    xdot    = delta_xdot
    xdotdot = delta_xdotdot
    y   = y0 + delta_y
    ydot    = delta_ydot
    ydotdot = delta_ydotdot
    z   = z0 + delta_z
    zdot    = delta_zdot
    zdotdot = delta_zdotdot
    # print('initial xyz: ', x,y,z)
    # print('initial delta xyz: ', delta_x,delta_y,delta_z)

    fx0 = m * thrust_x(x,zdot,xdotdot,omega,t)
    fy0 = m * thrust_y(y,ydotdot,omega,t)
    fz0 = m * thrust_z(z,zdotdot,xdot,omega,t)

    f_array   = np.array([[fx0,fy0,fz0]])
    X_array   = np.array([[x,y,z]])
    t_array   = np.array(t)
    mp_array  = np.array(mp)
    m_array = np.array(m0)
    mp_biprop_array = np.array(mp)

    #DeltaV maneuver 1
    # t += tburn
    thrust_deltaV1, mp_deltaV1 = vac_thrust(deltaV_A_0,Isp,m0,tburn,De=0,pe=0)
    m -= mp_deltaV1
    f_array = np.append(f_array,[[thrust_deltaV1,0,0]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp_deltaV1)
    m_array = np.append(m_array,m)
    mp_biprop_array = np.append(mp_biprop_array,mp_deltaV1)
    t_array  = np.append(t_array,t)

    #Proximity operations A:
    while x >= x0_A-1 and x < x1_A:
        Vx = Vx_A
        t += dt
        ta += dt
        delta_x_new, delta_xdot_new, delta_xdotdot_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
        delta_y_new, delta_ydot_new, delta_ydotdot_new = error_y(error_ym,omega,t)
        delta_z_new, delta_zdot_new, delta_zdotdot_new = error_z(error_xm,error_zm,error_zdot,omega,t)

        delta_x_rel       = delta_x_new - delta_x
        delta_xdot_rel    = delta_xdot_new - delta_xdot
        delta_xdotdot_rel = delta_xdotdot_new - delta_xdotdot

        delta_y_rel       = delta_y_new - delta_y
        delta_ydot_rel    = delta_ydot_new - delta_ydot
        delta_ydotdot_rel = delta_ydotdot_new - delta_ydotdot

        delta_z_rel       = delta_z_new - delta_z
        delta_zdot_rel    = delta_zdot_new - delta_zdot
        delta_zdotdot_rel = delta_zdotdot_new - delta_zdotdot

        x       = x0_A + Vx*ta + delta_x_rel
        xdot    = Vx + delta_xdot_rel
        xdotdot = delta_xdotdot_rel

        y       = delta_y_rel
        ydot    = delta_ydot_rel
        ydotdot = delta_ydotdot_rel

        z       = delta_z_rel
        zdot    = delta_zdot_rel
        zdotdot = delta_zdotdot_rel

        fx = m * thrust_x(x,zdot,xdotdot,omega,t) * margin
        fy = m * thrust_y(y,ydotdot,omega,t) * margin
        fz = m * thrust_z(z,zdotdot,xdot,omega,t) * margin
        ftot = abs(fx + fy + fz)
        mp = act.RCSpropellant(ftot,dt,Isp)
        m -= mp

        f_array  = np.append(f_array,[[fx,fy,fz]],axis=0)
        X_array  = np.append(X_array,[[x,y,z]],axis=0)
        mp_array = np.append(mp_array,mp)
        m_array = np.append(m_array,m)
        t_array  = np.append(t_array,t)

        delta_x       = delta_x_new
        delta_xdot    = delta_xdot_new
        delta_xdotdot = delta_xdotdot_new

        delta_y       = delta_y_new
        delta_ydot    = delta_ydot_new
        delta_ydotdot = delta_ydotdot_new

        delta_z       = delta_z_new
        delta_zdot    = delta_zdot_new
        delta_zdotdot = delta_zdotdot_new

    #Delta V maneuver 2
    # t += tburn
    thrust_deltaV2, mp_deltaV2 = vac_thrust(deltaV_A_1,Isp,m,tburn,De=0,pe=0)
    m -= mp_deltaV2
    f_array = np.append(f_array,[[thrust_deltaV2,0,0]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp_deltaV2)
    m_array = np.append(m_array,m)
    mp_biprop_array = np.append(mp_biprop_array,mp_deltaV2)
    t_array  = np.append(t_array,t)

    #Delta V maneuver 3
    # t += tburn
    thrust_deltaV3, mp_deltaV3 = vac_thrust(deltaV_B_0,Isp,m0,tburn,De=0,pe=0)
    m -= mp_deltaV3
    f_array = np.append(f_array,[[thrust_deltaV3,0,0]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp_deltaV3)
    m_array = np.append(m_array,m)
    mp_biprop_array = np.append(mp_biprop_array,mp_deltaV3)
    t_array  = np.append(t_array,t)

    #Proximity operations B
    while x >= x0_B-1 and x < x1_B:
        t  += dt
        tb += dt
        Vx = Vx_B

        delta_x_new, delta_xdot_new, delta_xdotdot_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
        delta_y_new, delta_ydot_new, delta_ydotdot_new = error_y(error_ym,omega,t)
        delta_z_new, delta_zdot_new, delta_zdotdot_new = error_z(error_xm,error_zm,error_zdot,omega,t)

        delta_x_rel       = delta_x_new - delta_x
        delta_xdot_rel    = delta_xdot_new - delta_xdot
        delta_xdotdot_rel = delta_xdotdot_new - delta_xdotdot

        delta_y_rel       = delta_y_new - delta_y
        delta_ydot_rel    = delta_ydot_new - delta_ydot
        delta_ydotdot_rel = delta_ydotdot_new - delta_ydotdot

        delta_z_rel       = delta_z_new - delta_z
        delta_zdot_rel    = delta_zdot_new - delta_zdot
        delta_zdotdot_rel = delta_zdotdot_new - delta_zdotdot


        x       = x0_B + Vx*tb + delta_x_rel
        xdot    = Vx + delta_xdot_rel
        xdotdot = delta_xdotdot_rel
        y       = delta_y_rel
        ydot    = delta_ydot_rel
        ydotdot = delta_ydotdot_rel

        z       = delta_z_rel
        zdot    = delta_zdot_rel
        zdotdot = delta_zdotdot_rel

        fx = m * thrust_x(x,zdot,xdotdot,omega,t) * margin
        fy = m * thrust_y(y,ydotdot,omega,t) * margin
        fz = m * thrust_z(z,zdotdot,xdot,omega,t) * margin
        ftot = abs(fx + fy + fz)
        mp = act.RCSpropellant(ftot,dt,Isp)
        m -= mp


        f_array = np.append(f_array,[[fx,fy,fz]],axis=0)
        X_array = np.append(X_array,[[x,y,z]],axis=0)
        mp_array = np.append(mp_array,mp)
        m_array = np.append(m_array,m)
        t_array = np.append(t_array,t)

        delta_x       = delta_x_new
        delta_xdot    = delta_xdot_new
        delta_xdotdot = delta_xdotdot_new

        delta_y       = delta_y_new
        delta_ydot    = delta_ydot_new
        delta_ydotdot = delta_ydotdot_new

        delta_z       = delta_z_new
        delta_zdot    = delta_zdot_new
        delta_zdotdot = delta_zdotdot_new

    #Delta V maneuver 4
    # t += tburn
    thrust_deltaV4, mp_deltaV4 = vac_thrust(deltaV_B_1,Isp,m0,tburn,De=0,pe=0)
    m -= mp_deltaV4
    f_array = np.append(f_array,[[thrust_deltaV4,0,0]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp_deltaV4)
    m_array = np.append(m_array,m)
    mp_biprop_array = np.append(mp_biprop_array,mp_deltaV4)
    t_array  = np.append(t_array,t)

    #Delta V maneuver 5
    # t += tburn
    thrust_deltaV5, mp_deltaV5 = vac_thrust(deltaV_d_0,Isp,m0,tburn,De=0,pe=0)
    m -= mp_deltaV5
    f_array = np.append(f_array,[[thrust_deltaV5,0,0]],axis=0)
    X_array  = np.append(X_array,[[x,y,z]],axis=0)
    mp_array = np.append(mp_array,mp_deltaV5)
    m_array = np.append(m_array,m)
    mp_biprop_array = np.append(mp_biprop_array,mp_deltaV5)
    t_array  = np.append(t_array,t)
    while x >= x0_d-10 and x < x1_d:
        Vx  = Vx_d
        td += dt
        t  += dt

        delta_x_new, delta_xdot_new, delta_xdotdot_new = error_x(error_xm,error_zm,error_xdot,error_zdot,omega,t)
        delta_y_new, delta_ydot_new, delta_ydotdot_new = error_y(error_ym,omega,t)
        delta_z_new, delta_zdot_new, delta_zdotdot_new = error_z(error_xm,error_zm,error_zdot,omega,t)

        delta_x_rel       = delta_x_new - delta_x
        delta_xdot_rel    = delta_xdot_new - delta_xdot
        delta_xdotdot_rel = delta_xdotdot_new - delta_xdotdot

        delta_y_rel       = delta_y_new - delta_y
        delta_ydot_rel    = delta_ydot_new - delta_ydot
        delta_ydotdot_rel = delta_ydotdot_new - delta_ydotdot

        delta_z_rel       = delta_z_new - delta_z
        delta_zdot_rel    = delta_zdot_new - delta_zdot
        delta_zdotdot_rel = delta_zdotdot_new - delta_zdotdot

        x       = x0_d + Vx*td + delta_x_rel
        xdot    = Vx + delta_xdot_rel
        xdotdot = delta_xdotdot_rel
        y       = delta_y_rel
        ydot    = delta_ydot_rel
        ydotdot = delta_ydotdot_rel

        z       = delta_z_rel
        zdot    = delta_zdot_rel
        zdotdot = delta_zdotdot_rel

        fx = m * thrust_x(x,zdot,xdotdot,omega,t) * margin
        fy = m * thrust_y(y,ydotdot,omega,t) * margin
        fz = m * thrust_z(z,zdotdot,xdot,omega,t) * margin
        ftot = abs(fx + fy + fz)
        mp = act.RCSpropellant(ftot,dt,Isp)
        m -= mp

        f_array = np.append(f_array,[[fx,fy,fz]],axis=0)
        X_array = np.append(X_array,[[x,y,z]],axis=0)
        mp_array = np.append(mp_array,mp)
        m_array = np.append(m_array,m)
        t_array = np.append(t_array,t)

        delta_x       = delta_x_new
        delta_xdot    = delta_xdot_new
        delta_xdotdot = delta_xdotdot_new

        delta_y       = delta_y_new
        delta_ydot    = delta_ydot_new
        delta_ydotdot = delta_ydotdot_new

        delta_z       = delta_z_new
        delta_zdot    = delta_zdot_new
        delta_zdotdot = delta_zdotdot_new



    #=================================================================================================================================================
    # Total thrust and propellant mass
    #=================================================================================================================================================
    #Delta V maneuvers in between phases
    thrust_deltaV = [thrust_deltaV1, thrust_deltaV2, thrust_deltaV3, thrust_deltaV4, thrust_deltaV5]
    mp_deltaV = mp_deltaV1 + mp_deltaV2 + mp_deltaV3 + mp_deltaV4 + mp_deltaV5

    #Total
    mp_tot     = np.sum(mp_array)

    #RCS thrust per engine
    RCS_thrust_x    = act.RCS_displacement_to_thrust(f_array[:,0],'y','normal')
    RCS_thrust_y    = act.RCS_displacement_to_thrust(f_array[:,1],'x','normal')
    RCS_thrust_z    = act.RCS_displacement_to_thrust(f_array[:,2],'z','normal')

    #=================================================================================================================================================
    # Errors
    #=================================================================================================================================================
    #========================================================================
    # Thruster misalignment
    #========================================================================
    #Misalignment
    angle = 2*np.pi/180

    T_error_y1, T_error_x1, T_error_z1 = act.thrust_error(max(f_array[:,0]),cg_orbit,angle)
    T_error_y2, T_error_x2, T_error_z2 = act.thrust_error(max(f_array[:,1]),cg_orbit,angle)
    T_error_y3, T_error_x3, T_error_z3 = act.thrust_error(max(f_array[:,2]),cg_orbit,angle)
    T_error_x = T_error_x1 + T_error_x2 + T_error_x3
    T_error_y = T_error_y1 + T_error_y2 + T_error_y3
    T_error_z = T_error_z1 + T_error_z2 + T_error_z3


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
    T_error_x_tot = T_error_x + Tgx + Tsp + Tm
    T_error_y_tot = T_error_y + Tgy
    T_error_z_tot = T_error_z

    RCS_thrust_error_x = margin * act.RCS_torque_to_thrust(T_error_y_tot,'y',cg_orbit,'error_bottom')
    RCS_thrust_error_y = margin * act.RCS_torque_to_thrust(T_error_x_tot,'x',cg_orbit,'error_bottom')
    RCS_thrust_error_z = margin * act.RCS_torque_to_thrust(T_error_z_tot,'z',cg_orbit,'error_bottom')

    RCS_thrust_error_tot = max(RCS_thrust_error_x,RCS_thrust_error_y,RCS_thrust_error_z)

    RCS_impulse_error = (np.sum(RCS_thrust_error_x*t) + np.sum(RCS_thrust_error_y*t) + np.sum(RCS_thrust_error_z*t))
    mp_error           = RCS_impulse_error / (Isp * g_earth)

    #=================================================================================================================================================
    # Rotations
    #=================================================================================================================================================


    angle_transfer      = 180 * np.pi / 180
    angle_rendezvous    = 180 * np.pi / 180

    time_transfer       = time_hohmann
    time_rendezvous     = time_hohmann
    RCS_thrust_rotations  = max(RCS_thrust_x)
    tburn                 = tburn_rot


    transfer_time,RCS_torque_rotations        = act.slew(RCS_thrust_rotations,tburn,angle_transfer,Iy)
    mp_rotations                 = 4 * 3 * act.RCSpropellant(RCS_thrust_rotations,tburn,Isp)

    T_error_rot_y, T_error_rot_x, T_error_rot_z = act.thrust_error(RCS_thrust_rotations,cg_orbit,angle)
    T_error_rot_x += Tgx + Tsp + Tm
    RCS_error_rot_x  = margin * act.RCS_torque_to_thrust(T_error_rot_x,'z',cg_orbit,'error_bottom')
    RCS_error_rot_y  = margin * act.RCS_torque_to_thrust(T_error_rot_y,'z',cg_orbit,'error_bottom')
    RCS_error_rot_z  = margin * act.RCS_torque_to_thrust(T_error_rot_z,'z',cg_orbit,'error_bottom')
    RCS_error_rot    = max(RCS_error_rot_x, RCS_error_rot_y, RCS_error_rot_z)


    #=================================================================================================================================================
    # Final thrust and propellant values
    #=================================================================================================================================================
    RCS_docking_min   = [min(RCS_thrust_x), min(RCS_thrust_y), min(RCS_thrust_z)]
    print(RCS_docking_min[0],RCS_docking_min[1],RCS_docking_min[2])
    F_docking_min     = [4*RCS_docking_min[0], 2*RCS_docking_min[1], 2*RCS_docking_min[2]]
    RCS_docking_max   = [max(RCS_thrust_x), max(RCS_thrust_y), max(RCS_thrust_z)]
    F_docking_max     = [4*RCS_docking_max[0], 2*RCS_docking_max[1], 2*RCS_docking_max[2]]
    mp_docking        = mp_tot
    t_docking         = t_array[-1]

    T_error_docking   = max(RCS_thrust_error_x,RCS_thrust_error_y,RCS_thrust_error_z)
    RCS_error_docking = max([RCS_thrust_error_x, RCS_thrust_error_y, RCS_thrust_error_z])
    mp_error_docking  = mp_error

    T_rot             = RCS_torque_rotations
    RCS_rot           = RCS_thrust_rotations
    mp_rot            = mp_rotations
    t_rot             = transfer_time

    T_error_rot       = [T_error_rot_x, T_error_rot_y, T_error_rot_z]
    RCS_error_rot     = max(RCS_error_rot_x, RCS_error_rot_y, RCS_error_rot_z)
    mp_error_rot      = act.RCSpropellant(RCS_error_rot,tburn,Isp)

    mp_total          = mp_docking + mp_error_docking + mp_rot + mp_error_rot
    acc_lat           = act.RCS_thrust_to_torque(F_docking_max[0],'y',cg_orbit)/Iy
    acc_long          = act.RCS_thrust_to_torque(F_docking_max[0],'z',cg_orbit)/Iz

    if printing == True:

        print('==========================================================')
        print('==========================================================')
        print('Docking time:',ta/60,tb/60,td/60,t_array[-1]/60, (transfer_time*2)/60, 2*period_phasing/60,t_array[-1]/60+ (transfer_time*2)/60+2*period_phasing/60)

        print('')
        print('Thrust for initial and final delta Vs: ', thrust_deltaV)
        print('')
        print('DOCKING')
        print('Max thrust total             : ', F_docking_max)
        print('Max thrust per engine        : ', RCS_docking_max)
        print('Min thrust total             : ', F_docking_min)
        print('Min thrust per engine        : ', RCS_docking_min)
        print('Propellant needed            : ', mp_docking)
        print('time                         : ', t_docking)
        print('REDUNDANCY')
        print('Disturbance torque           : ', T_error_docking)
        print('Redundancy thrust per engine : ', RCS_error_docking)
        print('Propellant needed            : ', mp_error_docking)
        print('')
        print('RV')
        print('Total torque                 : ', T_rot)
        print('Thrust                       : ', RCS_rot)
        print('Propellant needed            : ', mp_rot)
        print('REDUNDANCY')
        print('Disturbance torque           : ', T_error_rot)
        print('Redundancy thrust            : ', RCS_error_rot)
        print('Propellant used              : ', mp_error_rot)
        print('time                         : ', t_rot)
        print('time limit                   : ', time_rendezvous)
        print('')
        print('Total propellant needed      : ', mp_total)
        print('==========================================================')
        print('REQUIREMENTS: 0.16 m/s^2 acceleration')
        print('Acceleration longitudinal    :',acc_lat)
        print('Acceleration lateral         :',acc_long)
        print('==========================================================')
        print('==========================================================')
    #=====================================================================================================================================================================================================================
    # Plotting
    #=====================================================================================================================================================================================================================
    if plotting:
        fig, axs = plt.subplots(2, 3, constrained_layout=True)
        fig.suptitle('Relative position wrt target',fontsize=16)
        #Fx vs time
        axs[0][0].plot(t_array[2:],X_array[2:,0],color="navy")
        axs[0][0].grid(color="gainsboro")
        axs[0][0].set_xlabel("Time [s]")
        axs[0][0].set_ylabel("X position [m]")

        #Fx vs time
        axs[0][1].plot(t_array[2:],X_array[2:,1],color="navy")
        axs[0][1].grid(color="gainsboro")
        axs[0][1].set_xlabel("Time [s]")
        axs[0][1].set_ylabel("Y position [m]")

        #Fx vs time
        axs[0][2].plot(t_array[2:],X_array[2:,2],color="navy")
        axs[0][2].grid(color="gainsboro")
        axs[0][2].set_xlabel("Time [s]")
        axs[0][2].set_ylabel("Z postion [m]")
        plt.show()

        fig, axs = plt.subplots(2, 2, constrained_layout=True)
        fig.suptitle('RCS thrust',fontsize=16)
        #Fx vs time
        axs[0][0].plot(t_array,f_array[:,0],color="navy")
        axs[0][0].grid(color="gainsboro")
        axs[0][0].set_xlabel("Time [s]")
        axs[0][0].set_ylabel("Thrust in X [N]")

        #Fx vs time
        axs[0][1].plot(t_array,f_array[:,1],color="navy")
        axs[0][1].grid(color="gainsboro")
        axs[0][1].set_xlabel("Time [s]")
        axs[0][1].set_ylabel("Thrust in Y [N]")

        #Fx vs time
        axs[1][0].plot(t_array,f_array[:,2],color="navy")
        axs[1][0].grid(color="gainsboro")
        axs[1][0].set_xlabel("Time [s]")
        axs[1][0].set_ylabel("Thrust in Z [N]")

        #Propellant mass
        axs[1][1].plot(t_array,mp_array,color="navy")
        axs[1][1].grid(color="gainsboro")
        axs[1][1].set_xlabel("Time [s]")
        axs[1][1].set_ylabel("Propellant mass [kg]")
        plt.show()

        plt.figure()
        plt.plot(t_array,(m0-m_array),color="navy")
        plt.grid(color="gainsboro")
        plt.title("Time vs Propellant mass")
        plt.xlabel("Time [s]")
        plt.ylabel("Propellant mass [kg]")
        plt.show()

    return acc_lat, acc_long, mp_total, mp_rot+mp_error_rot, mp_docking+mp_error_docking, RCS_rot+RCS_error_rot,RCS_docking_max+RCS_error_docking,t_rot
