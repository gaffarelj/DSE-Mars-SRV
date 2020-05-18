from tradeoff import tradeoff_class as tc
from mass import Detailed_Mass_Budget as dmb
from complexity import ssto as c1
from complexity import ms as c2
from complexity import spaceplane as c3
from complexity import elevator as c4
import multiprocessing as mp
import numpy as np

dmass = tc.param(name="Dry Mass", weight=0.10963, direc="LB", Limitype="SD",Limit_val=1.5)
prmass = tc.param(name="Propellant Mass", weight=0.24875, direc="LB", Limitype="SD",Limit_val=1.5)
trl = tc.param(name="TRL", weight=0.0525)
comp = tc.param(name="Complexity", weight=0.28913, direc="LB", func="IRTS", Limitype="SD",Limit_val=1.5)
lov = tc.param(name="LOV Risk", weight=0.285, direc="LB", func="IRTS", Limitype="SD",Limit_val=1.2)
lom = tc.param(name="LOM Risk", weight=0.015, direc="LB", func="IRTS", Limitype="SD",Limit_val=1.2)

Isp = 400
d1_m = dmb.Mass_conc(5894.46, 0, Isp, "SSTO")
d1_in_list = [d1_m["dry"], d1_m["prop"], 5, c1.comp(), 1/38.3, 1/140.1152]
d1 = tc.design(name="Single Stage", sourcelist=d1_in_list)

d2_m = dmb.Mass_conc(1730.17, 3951.455, Isp, "2_stage")
d2_in_list = [d2_m[0], d2_m[1], 5, c2.comp(), 1/19.8203, 1/86.7670]
d2 = tc.design(name="Multi Stage", sourcelist=d2_in_list)

d3_m = dmb.Mass_conc(5023.57, 0, Isp, "SPACEPLANE")
d3_in_list = [d3_m["dry"], d3_m["prop"], 4, c3.comp(), 1/42.1071, 1/74.4518]
d3 = tc.design(name="Spaceplane", sourcelist=d3_in_list)

print(c4.comp())
d4_m = dmb.Mass_conc(0, 0, Isp, "SE")
d4_in_list = [d4_m[0], d4_m[1], 4, c4.comp(), 1/82.3690, 1/280.0512]
d4 = tc.design(name="Space Elevator", sourcelist=d4_in_list)

tradeoff = tc.tradeoff(design_list = [d1, d2, d3], param_list= [dmass, prmass, trl, comp, lov, lom])

tradeoff.get_tradeoff()
colors = [tc.color("EF5350", "red"), tc.color("FB8C00", "orange"), tc.color("FFEB3B", "yellow"), tc.color("8BC34A", "green"), tc.color("00BCD4", "blue")]
tradeoff.get_output(language="python", color_list=colors, width=13)
input()
sens = tc.sensitivity(tradeoff, samples=1)
#sens.addto_technical(0.25)
sens.addto_weights(0.5)
sens.get_RMS()
#sens.get_sens_linux()
#print(sens.per)

do_analysis = True
if __name__ == "__main__" and do_analysis:
	n_d = len(sens.tro.design_list)
	deltas = [[] for i in range(n_d)]
	n = int(5e4)
	ori = np.array([w.weight for w in sens.tro.param_list])
	for i in range(n):
		print(round(i/n*100, 2), end="\r")
		ret, wei = sens.sens(sens.n)
		wei = np.array(wei)
		delta = wei-ori
		deltas[np.where(ret==1)[0][0]].append(delta)
	print()
	weight_list = np.array([param.weight for param in sens.tro.param_list])
	print([p.name for p in sens.tro.param_list])
	for i in range(n_d):
		if len(deltas[i]) == 0:
			print("Concept never winner")
		else:
			print(np.average(deltas[i], axis=0))
			#print(np.average(deltas[i], axis=0)+2*np.std(deltas[i], axis=0))
			#print(np.average(deltas[i], axis=0)-2*np.std(deltas[i], axis=0))
