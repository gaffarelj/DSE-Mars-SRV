import tradeoff_class as tc

mass = tc.param(name="mass",weight=0.75)
param2 = tc.param(name="param2",weight=0.25)
design1 = tc.design(name="1",sourcelist=[1,2])
design2 = tc.design(name="2",sourcelist=[3,9])
design3 = tc.design(name="3",sourcelist=[4,6])
color1 = tc.color("FFFC9E","yellow")
color2 = tc.color("FE0000","red")



tradeoff =tc.tradeoff(design_list = [design1,design2,design3],param_list= [mass,param2],out="latex",color_list=[color1,color2])