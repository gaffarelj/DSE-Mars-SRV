import calc_class as cs

mech = cs.line(0.5,"bi")
fluid =  cs.line(1,"uni")
elec = cs. line(1,"uni")
data = cs.line(1,"bi")

valve = cs.system(1)
controller = cs.system(5)
pump = cs.system(2)
fil = cs.system(1)
motor = cs.system(3)

valve.add_c(mech,[controller, fil, pump])
valve.add_c(fluid,[fil, pump])
valve.add_c(elec,controller)
valve.add_c(data,controller)
controller.add_c(mech,motor)
controller.add_c(elec,motor)
pump.add_c(mech,motor)
pump.add_c(elec,motor)

result = cs.complexity()
print(result.structural, result.structural == 45.57519182318174)