import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import disturbances as dist
import actuator_properties as act
import reentry_control as rty
import rendez_vous as rvs
import ascent as asc
import pitching_moment as pm
import DYNAS_new as dns

import numpy as np
from matplotlib import pyplot as plt

Isp = act.Isp_mono

#=======================================================================================
#Reentry
#=======================================================================================
I_reentry         = [4053828.763, 372883.4811, 4053828.763]
cg_reentry        = act.z_cg_empty
tburn_reentry     = 1.
reentry_RCS_rot   = 451.85155228408263
pitching_time     = pm.time
pitching_moment   = pm.pitching_moment

reentry_mp, reentry_mp_rot, reentry_mp_pitch, reentry_RCS_rot, reentry_RCS_pitch, reentry_t_rot = rty.reentry_control(pitching_time,pitching_moment,I_reentry,cg_reentry, tburn_reentry, reentry_RCS_rot, Isp)


#=======================================================================================
#Ascent
#=======================================================================================
I_ascent         = [4321654.768, 407613.6669, 4321654.768]
cg_ascent        = act.z_cg_end_ascent
tburn_ascent     = 1.
ascent_RCS_rot   = 451.85155228408263
pitching_moment  = dns.am
pitching_time    = dns.t

ascent_mp, ascent_mp_rot, ascent_mp_pitch, ascent_RCS_rot, ascent_RCS_pitch, ascent_t_rot =asc.ascent_control(pitching_time,pitching_moment,I_ascent,cg_ascent, tburn_ascent, ascent_RCS_rot, Isp)

#=======================================================================================
#RV and Docking
#=======================================================================================
I_RV             = [4124615.776, 433619.6909, 4124615.776]
cg_RV            = act.z_cg_orbit
tburn_RV_deltaV  = 6.
tburn_RV_rot     = 1.
m0               = 53248.16053543461

#Navigation measurement errors:
error = [1.,1.,1.,1.,1.,1.]

RV_acc_lat, RV_acc_long, RV_mp, RV_mp_rot, RV_mp_docking, RV_RCS_rot, RV_RCS_docking, RV_t_rot = rvs.RV_Docking_control(I_RV,cg_RV,tburn_RV_deltaV,tburn_RV_rot,Isp,m0,error)


# rho_propellant = 1440.
# Mpropellant_total = reentry_mp + ascent_mp + RV_mp
# Vpropellant_total = Mpropellant_total / rho_propellant
# print(Mpropellant_total,Vpropellant_total)
