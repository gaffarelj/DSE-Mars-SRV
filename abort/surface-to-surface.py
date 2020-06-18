import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from astrodynamics import astro_tools_nonfull as AT
from astrodynamics import BAT


target = 3020
t_b = 120
tilt = 4
return_mass_decrease = 51650+6000
m_prop = 168953-return_mass_decrease
m = 36025+m_prop

mars=AT.Planet()
long0, lat0 = np.radians(-27.088), np.radians(4.51)

ascent = BAT.ascent_sim(tb=t_b,initial_tilt=tilt,i_base=41,h0=-3*10**3,d=6,M_initial=m,Mp_class2=m_prop,Isp=383.25,n=9,De=1.35,pe=6077.91)[4]
downr = np.sqrt(((ascent["long"]-long0)*mars.r)**2 + ((ascent["lat"]-lat0)*mars.r)**2)
print("downrange:", downr[-1]/1000, "[km]")
print("propellant used:", ascent["Mprop"][-1][0], "[kg]")
AT.plot_dual(ascent["time"], [h/1000 for h in ascent["altitude"]], [target-d/1000 for d in downr], "Time [s]", "Altitude [km]", "Target distance [km]")