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
d1 = tc.design(name="Singe stage", sourcelist=d1_in_list)

d2_m = dmb.Mass_conc(1730.17, 3951.455, Isp, "2_stage")
d2_in_list = [d2_m[0], d2_m[1], 5, c2.comp(), 1/19.8203, 1/86.7670]
d2 = tc.design(name="Multi stage", sourcelist=d2_in_list)

d3_m = dmb.Mass_conc(5023.57, 0, Isp, "SPACEPLANE")
d3_in_list = [d3_m["dry"], d3_m["prop"], 4, c3.comp(), 1/42.1071, 1/74.4518]
d3 = tc.design(name="Spaceplane", sourcelist=d3_in_list)


test = tc.tradeoff(design_list = [d1, d2, d3], param_list= [dmass, prmass, trl, comp, lov, lom])

test.get_tradeoff()


Rel = tc.param(name="Reliability", weight=0.30, Limitype="fixed",Limit_val=[0,1])
cost = tc.param(name="Cost", weight=0.35, Limitype="fixed",Limit_val=[0,1])
ruse = tc.param(name="Re-Usability", weight=0.20, Limitype="fixed",Limit_val=[0,1] )
sus = tc.param(name="Sustainability", weight=0.15, Limitype="fixed",Limit_val=[0,1])


o_list = np.ndarray([len(test.param_list),len(test.design_list)])
for i in range(len(test.param_list)):
    param = test.param_list[i]
    print(param.val_out)
    o_list[i] = np.array(param.val_out)
print(o_list.transpose())
w1 = np.array([0,0,0,0,0.95,0.05])
#w2 = np.array([0,0,0,0,0,0])
w2 = np.array([0.1275,0.125,0.5975,0.15,0,0])
w3 = np.array([0.1,0.50,0.4,0,0,0])
#w3 = np.array([0,0,0,0,0,0])
w4 = np.array([0.3,0.7,0,0,0,0])
#w4 = np.array([0,0,0,0,0,0])

p1 = np.dot(o_list.transpose(),w1)
print(p1)
p2 = np.dot(o_list.transpose(),w2)
p3 = np.dot(o_list.transpose(),w3)
p4 = np.dot(o_list.transpose(),w4)



d1 = tc.design(name="Singe stage", sourcelist=[p1[0],p2[0],p3[0],p4[0]])


d2 = tc.design(name="Multi stage", sourcelist=[p1[1],p2[1],p3[1],p4[1]])


d3 = tc.design(name="Spaceplane", sourcelist=[p1[2],p2[2],p3[2],p4[2]])


tnew = tc.tradeoff(design_list = [d1, d2, d3], param_list= [Rel, cost, ruse, sus])

tnew.get_tradeoff()
colors = [tc.color("EF5350", "red"), tc.color("FB8C00", "orange"), tc.color("FFEB3B", "yellow"), tc.color("8BC34A", "green"), tc.color("00BCD4", "blue")]
tnew.get_output()
#tnew.get_output(language="latex", color_list=colors, width=15)
