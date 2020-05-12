from tradeoff import tradeoff_class as tc
import Detailed_Mass_Budget as dmb
from complexity import ssto as c1
from complexity import ms as c2
from complexity import spaceplane as c3
from complexity import elevator as c4

dmass = tc.param(name="Dry Mass", weight=0.10963, direc="LB", Limitype="minmax")
prmass = tc.param(name="Propellant Mass", weight=0.24875, direc="LB", Limitype="minmax")
trl = tc.param(name="TRL", weight=0.0525)
comp = tc.param(name="Complexity", weight=0.28913, direc="LB", func="DRTS", Limitype="minmax")
lov = tc.param(name="LOV Risk", weight=0.285, direc="LB", func="IRTS", Limitype="minmax")
lom = tc.param(name="LOM Risk", weight=0.015, direc="LB", func="IRTS", Limitype="minmax")


d1_in_list = [*dmb.Mass_conc(5894.46, 0, 410, "SSTO"), 5, c1.comp(), 1/38.3, 1/140.1152]
d1 = tc.design(name="Singe stage", sourcelist=d1_in_list)

d2_in_list = [*dmb.Mass_conc(1730.17, 3951.455, 410, "2_stage"), 5, c2.comp(), 1/19.8203, 1/86.7670]
d2 = tc.design(name="Multi stage", sourcelist=d2_in_list)

d3_in_list = [*dmb.Mass_conc(5023.57, 0, 410, "SPACEPLANE"), 4, c3.comp(), 1/42.1071, 1/74.4518]
d3 = tc.design(name="Spaceplane", sourcelist=d3_in_list)

d4_in_list = [*dmb.Mass_conc(0, 0, 410, "SE"), 4, c4.comp(), 1/82.3690, 1/280.0512]
d4 = tc.design(name="Space Elevator", sourcelist=d4_in_list)


#if not ([d1_in_list, d2_in_list, d3_in_list, d4_in_list] == [[32251.860818340043, 149696.55314834046, 5, 5.185715438377493, 0.026109660574412535, 0.0071369844242451935], [23919.610356608104, 65190.20395038708, 5, 5.330214381772641, 0.050453323108126516, 0.011525118996853644], [28897.978921347618, 106800.0508397109, 4, 5.524977891358354, 0.023748963951447617, 0.013431508707647094], [3853913.6882668925, 179427.97949182583, 4, 5.963027375592036, 0.012140489747356408, 0.0035707756295991593]]):
#	raise Exception("Mass budget changed !")


tradeoff =tc.tradeoff(design_list = [d1, d2, d3], param_list= [dmass, prmass, trl, comp, lov, lom])

tradeoff.get_tradeoff()
col1 = tc.color("EF5350", "red")
col2 = tc.color("FB8C00", "orange")
col3 = tc.color("FFEB3B", "yellow")
col4 = tc.color("8BC34A", "green")
col5 = tc.color("00BCD4", "blue")
tradeoff.get_output(language="latex", color_list=[col1, col2, col3, col4, col5], width=15)
"""
sens = tc.sensitivity(tradeoff, 10000)
sens.addto_technical(0.1)
sens.addto_weights(0.1)
sens.get_RMS()
sens.get_sens()
print(sens.per)
print(sens.RMS)"""