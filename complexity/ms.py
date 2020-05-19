from complexity import calc_class as cs
def comp():
    cs.reset()
    power = cs.line(0.25,"uni")
    data = cs.line(1,"bi")
    structural = cs.line(0.1,"bi")
    flow = cs.line(1,"uni")
    hydra = cs.line(0.3,"bi")

    locking_mechanism = cs.system(1)
    hatch = cs.system(1)
    deceleration_system = cs.system(3)
    s2_gn_comp = cs.system(2)
    s2_sensors = cs.system(1.5)
    bio_monitoring = cs.system(0.5)
    f_comp = cs.system(2)
    enc_dec = cs.system(0.5)
    amplifiers = cs.system(0.5)
    receving_ant = cs.system(0.5)
    transmitting = cs.system(0.5)
    s2_batteries = cs.system(1.5)
    fuel_cells = cs.system(3.5)
    s2_pcdu = cs.system(0.5)
    heater = cs.system(1)
    atmos = cs.system(3)
    oxygentanks = cs.system(0)
    cargo_compartment = cs.system(1)
    heat_shield = cs.system(4)
    cooling_system = cs.system(3)
    s2_actuator_tanks = cs.system(0)
    s2_actuator_feedsystem = cs.system(1.5)
    s2_control_actuators = cs.system(2)
    s2_tanks = cs.system(1)
    s2_hydraulics = cs.system(1.5)
    s2_chamber= cs.system(2.5)
    s2_feedsystem = cs.system(3.5)
    s2_pressurization = cs.system(3.5)
    s2_landing_legs = cs.system(1)
    s1_struc = cs.system(4.5)
    s2_struc = cs.system(4)

    s1_batteries = cs.system(1.5)
    s1_hydraulics = cs.system(1.5)
    s1_chamber= cs.system(2.5)
    s1_feedsystem = cs.system(3.5)
    s1_pressurization = cs.system(3.5)
    s1_landing_legs = cs.system(1)
    s1_actuator_tanks = cs.system(0)
    s1_actuator_feedsystem = cs.system(1.5)
    s1_control_actuators = cs.system(2)
    s1_pcdu = cs.system(0.5)
    s1_gn_comp = cs.system(2)
    s1_sensors = cs.system(1.5)
    s1_tanks = cs.system(1)

    stage = cs.system(2)

    locking_mechanism.add_c(structural,s2_struc)
    hatch.add_c(structural,s2_struc)
    deceleration_system.add_c(structural,s2_struc)
    s2_batteries.add_c(structural,s2_struc)
    s2_batteries.add_c(power,s2_pcdu)
    fuel_cells.add_c(structural,s2_struc)
    fuel_cells.add_c(power,s2_pcdu)
    oxygentanks.add_c(structural,s2_struc)
    oxygentanks.add_c(flow,atmos)
    atmos.add_c(data,heater)
    s2_sensors.add_c(data,s2_gn_comp)
    s2_gn_comp.add_c(data,[deceleration_system, f_comp, s2_landing_legs, s2_feedsystem, s2_actuator_feedsystem])
    s2_pcdu.add_c(power,[heater, atmos, amplifiers, s2_gn_comp, f_comp])
    bio_monitoring.add_c(data,[f_comp, enc_dec, atmos])
    f_comp.add_c(data,[cooling_system, enc_dec])
    enc_dec.add_c(data,amplifiers)
    amplifiers.add_c(data,[receving_ant, transmitting])
    cargo_compartment.add_c(structural,s2_struc)
    heat_shield.add_c(structural,s2_struc)
    cooling_system.add_c(flow,heat_shield)
    s2_actuator_tanks.add_c(structural,s2_struc)
    s2_actuator_tanks.add_c(flow,s2_actuator_feedsystem)
    s2_control_actuators.add_c(structural,s2_struc)
    s2_actuator_feedsystem.add_c(flow,s2_control_actuators)
    s2_tanks.add_c(structural,s2_struc)
    s2_tanks.add_c(flow,s2_pressurization)
    s2_hydraulics.add_c(structural,s2_struc)
    s2_hydraulics.add_c(hydra,[s2_control_actuators, s2_feedsystem, s2_chamber, s2_landing_legs])
    s2_pressurization.add_c(flow,s2_feedsystem)
    s2_feedsystem.add_c(flow,s2_chamber)
    s2_chamber.add_c(structural,s2_struc)
    s2_landing_legs.add_c(structural,s2_struc)

    stage.add_c(structural,[s1_struc, s2_struc])
    stage.add_c(data,s2_gn_comp)

    s2_gn_comp.add_c(data,s1_gn_comp)

    s1_batteries.add_c(structural,s1_struc)
    s1_batteries.add_c(power,s1_pcdu)
    s1_sensors.add_c(data,s1_gn_comp)
    s1_gn_comp.add_c(data,[s1_landing_legs, s1_feedsystem])
    s1_gn_comp.add_c(power,s1_pcdu)
    s1_gn_comp.add_c(data,s1_actuator_feedsystem)
    s1_actuator_tanks.add_c(structural,s1_struc)
    s1_actuator_tanks.add_c(flow,s1_actuator_feedsystem)
    s1_control_actuators.add_c(structural,s1_struc)
    s1_actuator_feedsystem.add_c(flow,s1_control_actuators)
    s1_hydraulics.add_c(structural,s1_struc)
    s1_hydraulics.add_c(hydra,[s1_control_actuators, s1_feedsystem, s1_chamber, s1_landing_legs])
    s1_pressurization.add_c(flow,s1_feedsystem)
    s1_feedsystem.add_c(flow,s1_chamber)
    s1_chamber.add_c(structural,s1_struc)
    s1_landing_legs.add_c(structural,s1_struc)
    s1_tanks.add_c(structural,s1_struc)
    s1_tanks.add_c(flow,s1_pressurization)

    result = cs.complexity()
    return result.average
    print(result.structural, result.structural == 244.54678488250926)
    print(result.average, result.average == 5.434372997389095)