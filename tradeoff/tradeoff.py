import tradeoff_class as tc
def f1():
    return 1
def f2():
    return 2
def f3():
    return 3
def f4():
    return 4
def f5():
    return 5
def f6():
    return 6


mass = tc.param(name="mass",weight=0.75)
param2 = tc.param(name="param2",weight=0.25)
design1 = tc.design(name="1",sourcelist=[f1,f4])
design2 = tc.design(name="2",sourcelist=[f2,f5])
design3 = tc.design(name="3",sourcelist=[f3,f6])
color1 = tc.color("FFFC9E","yellow")
color2 = tc.color("FE0000","red")



tradeoff =tc.tradeoff(design_list = [design1,design2,design3],param_list= [mass,param2],out="latex",color_list=[color1,color2])