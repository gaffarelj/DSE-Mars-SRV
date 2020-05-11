import tradeoff_class as tc

mass = tc.param(name="mass",weight=0.75,Limitype="2SD")
param2 = tc.param(name="param2",weight=0.25,Limitype="2SD")
design1 = tc.design(name="1",sourcelist=[1,20])
design2 = tc.design(name="2",sourcelist=[3,90])
design3 = tc.design(name="3",sourcelist=[4,60])
color1 = tc.color("FFFC9E","yellow")
color2 = tc.color("FE0000","red")



tradeoff =tc.tradeoff(design_list = [design1,design2,design3],param_list= [mass,param2])

tradeoff.get_tradeoff()
tradeoff.get_output(language="python",color_list=[color1,color2])

sens = tc.sensitivity(tradeoff,20000)
sens.addto_technical(0.0025)
sens.get_sens()
print(sens.per)