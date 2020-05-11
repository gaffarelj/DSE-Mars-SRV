from tradeoff import tradeoff_class as tc
import Detailed_Mass_Budget as dmb
from complexity import ssto as c1
from complexity import ms as c2
from complexity import spaceplane as c3
from complexity import elevator as c4

dmass = tc.param(name="Dry Mass",weight=0.10963,direc="LB")
prmass = tc.param(name="Propellant Mass",weight=0.24875,direc="LB")
trl = tc.param(name="TRL",weight=0.0525)
comp = tc.param(name="Complexity",weight=0.28913,direc="LB",func="IRTS")
lov = tc.param(name="LOV Risk",weight=0.285,direc="LB",func="IRTS")
lom = tc.param(name="LOM Risk",weight=0.015,direc="LB",func="IRTS")


d1_in_list = [dmb.Mass_conc(5894.46,0,410,"SSTO")[0][0],dmb.Mass_conc(5894.46,0,410,"SSTO")[0][1],5,c1.comp(),1/38.3,1/140.1152]
d1 = tc.design(name="Singe stage",sourcelist=d1_in_list)

d2_in_list = [dmb.Mass_conc(1730.17,3951.455,410,"2_stage")[0][0],dmb.Mass_conc(730.17,3951.455,410,"2_stage")[0][1],5,c2.comp(),1/19.8203,1/86.7670]
d2 = tc.design(name="Multi stage",sourcelist=d2_in_list)

d3_in_list = [dmb.Mass_conc(5023.57,0,410,"SPACEPLANE")[0][0],dmb.Mass_conc(5023.57,0,410,"SPACEPLANE")[0][1],4,c3.comp(),1/42.1071,1/74.4518]
d3 = tc.design(name="Spaceplane",sourcelist=d3_in_list)

d4_in_list = [dmb.Mass_conc(0,0,410,"SE")[0][0],dmb.Mass_conc(0,0,410,"SE")[0][1],4,c4.comp(),1/82.3690,1/280.0512]
d4 = tc.design(name="Space Elevator",sourcelist=d4_in_list)

color1 = tc.color("FFFC9E","yellow")
color2 = tc.color("FE0000","red")



tradeoff =tc.tradeoff(design_list = [d1,d2,d3],param_list= [dmass,prmass,trl,comp,lov,lom])

tradeoff.get_tradeoff()
tradeoff.get_output()
#tradeoff.get_output(language="latex",color_list=[color1,color2])
sens = tc.sensitivity(tradeoff,10000)
sens.addto_technical(0.1)
sens.addto_weights(0.1)
sens.get_RMS()
sens.get_sens()
print(sens.per)
print(sens.RMS)