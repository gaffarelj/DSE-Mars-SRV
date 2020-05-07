import calc_class as cs

power = cs.line(0.25,"uni")
data = cs.line(1,"bi")
structural = cs.line(0.1,"bi")
structural_spec = cs.line(1.5,"bi")
flow = cs.line(1,"uni")
laser = cs.line(2,"bi")
hydra = cs.line(0.3,"bi")

locking_mechanism = cs.system(1)
hatch = cs.system(1)
deceleration_system = cs.system(3)
gn_comp = cs.system(2)
sensors = cs.system(1.5)
launch_abort = cs.system(1.5)
bio_monitoring = cs.system(0.5)
f_comp = cs.system(2)
enc_dec = cs.system(0.5)
amplifiers = cs.system(0.5)
receving_ant = cs.system(0.5)
transmitting = cs.system(0.5)
pcdu = cs.system(0.5)
heater = cs.system(1)
atmos = cs.system(3)
oxygentanks = cs.system(0)
cargo_compartment = cs.system(1)
cooling_system = cs.system(4)
capulse_struc = cs.system(5)
tether = cs.system(5)
anchor_point = cs.system(5)
counterweight = cs.system(5)
mirror = cs.system(4)
laser_cooling = cs.system(5)
Laser_generator = cs.system(5)
power_delivery = cs.system(3)
energy_reciever = cs.system(4)


l = [locking_mechanism,
hatch,
deceleration_system,
gn_comp,
sensors,
launch_abort,
bio_monitoring,
f_comp,
enc_dec,
amplifiers,
receving_ant,
transmitting,
pcdu,
heater,
atmos,
oxygentanks,
cargo_compartment,
cooling_system,
capulse_struc,
tether,
anchor_point,
counterweight,
mirror,
laser_cooling,
Laser_generator,
power_delivery,
energy_reciever]


locking_mechanism.add_c(structural,capulse_struc)
hatch.add_c(structural,capulse_struc)
deceleration_system.add_c(structural,capulse_struc)
launch_abort.add_c(structural,capulse_struc)
oxygentanks.add_c(structural,capulse_struc)
oxygentanks.add_c(flow,atmos)
atmos.add_c(data,heater)
sensors.add_c(data,gn_comp)
gn_comp.add_c(data,deceleration_system)
gn_comp.add_c(data,launch_abort)
gn_comp.add_c(data,f_comp)
pcdu.add_c(power,heater)
pcdu.add_c(power,atmos)
pcdu.add_c(power,amplifiers)
pcdu.add_c(power,gn_comp)
pcdu.add_c(power,f_comp)
bio_monitoring.add_c(data,f_comp)
bio_monitoring.add_c(data,enc_dec)
bio_monitoring.add_c(data,atmos)
f_comp.add_c(data,cooling_system)
f_comp.add_c(data,enc_dec)
enc_dec.add_c(data,amplifiers)
amplifiers.add_c(data,receving_ant)
amplifiers.add_c(data,transmitting)
cooling_system.add_c(flow,mirror)
anchor_point.add_c(structural_spec,tether)
tether.add_c(structural_spec,counterweight)
tether.add_c(structural_spec,capulse_struc)
energy_reciever.add_c(structural,capulse_struc)
laser_cooling.add_c(hydra,Laser_generator)
power_delivery.add_c(power,Laser_generator)
Laser_generator.add_c(laser,energy_reciever)
mirror.add_c(structural,capulse_struc)
energy_reciever.add_c(power,pcdu)

result = cs.complexity(l)
print(result.structural)
print(result.average)