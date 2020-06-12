import sys
sys.path.append('../astrodynamics')
import mars_standard_atmosphere as MSA
import DYNAS as dns
import disturbances as dist
import actuator_properties as act
import reentry_control as rty
import rendez_vous as rvs
import ascent as asc

import numpy as np
from matplotlib import pyplot as plt



rho_propellant = 1440.

Mpropellant_total = rvs.mp_total + asc.mp_total + rty.mp
Vpropellant_total = Mpropellant_total / rho_propellant


print(Mpropellant_total,Vpropellant_total)
