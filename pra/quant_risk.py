import quant_class as qc
import numpy as np
# Enter your events here as follows
# Basic failure events of a component as such
# "<id>" : qc.event("<description>",<failure probability for on component>,"<consequence from list>",count=<number of components>,redundancy=<number of redudant components>)
# Failure event that require a combination of events as
# "<id>" : qc.comp("<description>","<consequence from list>",<list of id's contributing to event>)
event_dict ={
    "e1_s" : qc.event("Propellant tank rupture",8.2e-6,"LOC",count= 2, source="fta:launchvehicle_risk"),
    "e2_s" : qc.event("Thrust structure failure",3.4e-6,"LOV",count= 9, redundancy=2, source="fta:launchvehicle_risk"),
    "e3_s" : qc.event("Pipe leakage",1.1e-5,"LOM",count= 1, source="fta:leakage_probability"),
	"e4_s" : qc.event("Landing leg collapsing",5e-5,"LOV",count=8, redundancy=1, source="fta:touchdown"),
    "e5_s" : qc.event("Skirt buckling",3.4e-6,"LOM",count=2, source="fta:launchvehicle_risk"),
    "e6_s" : qc.event("Separation mechanism failure",5e-4,"LOM",count= 1, source="fta:launchvehicle_risk"),
    "e8_s" : qc.event("Hydraulics system failure",6.96e-5,"LOV",count= 9, redundancy=2, source="fta:hydrau-relia"),
    "e9_s" : qc.event("Legs not deploying",1.67e-3,"LOV",count=8, redundancy=1, source="fta:spacex-leg"),
	"e10_s" : qc.event("Depressurisation of the capsule",7.23e-5,"LOV",count= 1, source="fta:orbit-debris"),
    "e11_s" : qc.event("Buckling of the capsule",3.40e-6,"LOM",count= 1, source="fta:launchvehicle_risk"),
    "e12_s" : qc.event("Attenuation system failure",6.96e-5,"LOV",count= 1, source="fta:hydrau-relia"),
	"e13_s" : qc.event("Side Hatch not opening",0,"LOM",count= 1, source="~Negligible"),
    "e14_s" : qc.event("Docking hatch not opening",0,"LOM",count= 1, source="~Negligible"),
	"e15_s" : qc.event("Capsule ring not deploying",6.96e-5,"LOM",count= 1, source="fta:hydrau-relia"),
    "e16_s" : qc.event("Not sealed connection",0,"LOM",count= 1, source="~Negligible"),

    "e1_th" : qc.event("TPS damaged due to surface debris",3.00e-6,"LOM",count= 1, source="fta:tps-fail"),
    "e2_th" : qc.event("TPS damaged due to space debris",7.23e-4,"LOC",count= 1, source="fta:orbit-debris"),
    "e3_th" : qc.event("TPS breaks due to fatigue",1.20e-7,"LOC",count= 1, source="fta:tps-fail"),
	"e4_th" : qc.event("TPS breaks due to high stress",1.44e-8,"LOC",count= 1, source="fta:tps-fail"),

	"e1_en" : qc.event("Catastrophic engine failure",1.09e-3,"LOV",count= 9, redundancy=2, source="fta:launchvehicle_risk"),
    "e2_en" : qc.event("Control system failure",1.37e-4,"LOV",count= 9, redundancy=2, source="fta:ren_fail"),
    "e3_en" : qc.event("Ignition failure",2.27e-5,"LOV",count= 9, redundancy=2, source="fta:ren_fail"),
	"e4_en" : qc.event("Fluid leak",1.10e-5,"LOM",count= 9, redundancy=2, source="fta:leakage_probability"),
    "e5_en" : qc.event("Fuel pump failure",3.33e-6,"LOM",count= 9, redundancy=2, source="fta:ESREL"),
    "e6_en" : qc.event("Oxidise pump failure",3.33e-6,"LOM",count= 9, redundancy=2, source="fta:ESREL"),
    "e7_en" : qc.event("Turbine failure",1.34e-6,"LOM",count= 9, redundancy=2, source="fta:gas_turbine_reliability"),
    "e8_en" : qc.event("Gimbal bearing failure",5.50e-4,"LOV",count= 9, redundancy=2, source="fta:ren_fail"),
    "e9_en" : qc.event("Combustion Instability",3.03e-3,"LOV",count= 9, redundancy=2, source="fta:ren_fail"),

	"e1_gnc" : qc.event("Docking camera failure",8.04e-5,"LOM",count = 4,redundancy = 2, source="~\\protect\\footnotemark[1] \\footnotetext[1]{Blackfly S Datasheet, \\url{https://flir.app.boxcn.net/s/mj27am7zik371ivyzv352gmt390yqt6z}, accessed on the 16th of June 2020.}"),
    "e2_gnc" : qc.event("Landing camera failure",8.04e-5,"LOM",count = 4,redundancy = 2, source="~\\protect\\footnotemark[1]"),
    "e3_gnc" : qc.event("Star tracker failure",7.5e-6,"LOV",count = 3,redundancy = 2, source="~\\protect\\footnotemark[1]"),
    "e4_gnc" : qc.event("Gyroscope failure",1.29e-3,"LOV",count = 4,redundancy = 3, source="fta:mtbf-risk"),
    "e5_gnc" : qc.event("Accelerometer failure",1.29e-3,"LOC",count = 4,redundancy = 3, source="fta:mtbf-risk"),
    "e6_gnc" : qc.event("Flight computer failure",7.41e-4,"LOC",count = 4,redundancy = 3, source="fta:launchvehicle_risk"),
    "e7_gnc" : qc.event("RCS Thruster failure",4.00e-3,"LOV",count = 32,redundancy = 6, source="fta:launchvehicle_risk"),
    "e8_gnc" : qc.event("Ignition failure",2.27e-5,"LOM",count = 1, source="fta:ren_fail"),
    "e9_gnc" : qc.event("Valve/Feedsystem Failure",5.60e-4,"LOV",count = 32,redundancy = 6, source="fta:utm-reliability"),
    "e10_gnc" : qc.event("Critical thruster failure:",1.09e-3,"LOV",count = 32,redundancy = 2, source="fta:launchvehicle_risk"),
    "e11_gnc" : qc.event("Propellant leakage",1.10e-5,"LOV",count = 2, source="fta:leakage_probability"),
    "e12_gnc" : qc.event("Control algorithm errors",3.71e-5,"LOM",count = 1, source="fta:launchvehicle_risk"),
    "e13_gnc" : qc.event("Guidance algorithm errors",3.71e-5,"LOM",count = 1, source="fta:launchvehicle_risk"),
    "e14_gnc" : qc.event("Engine thruster misalignment",0,"LOM",count = 32, source="~Negligible"),

    "e2_com" : qc.event("Antenna failure",8.69e-3,"LOM",count = 4,redundancy = 2, source="fta:reliability_communication"),
    "e3_com" : qc.event("Transmitter/receiver failure",2.53e-2,"LOM",count = 2, redundancy=1, source="fta:reliability_communication"),
    "e5_com" : qc.event("No line of sight post abort landing",5.00e-3,"LOM",count = 1, source="~Estimation"),

	"e1_pow" : qc.event("PDU failure",6.25e-4,"LOV",count = 4, redundancy=3, source="fta:failure-power"),
    "e2_pow" : qc.event("Fuel cell failure",2.50e-4,"LOM",count = 4,redundancy = 3, source="fta:failure-power"),
    "e3_pow" : qc.event("ICE overheating",1.25e-3,"LOM",count = 1, source="fta:failure-power"),
    "e4_pow" : qc.event("ICE failure",1.50e-3,"LOM",count = 2, redundancy=1, source="fta:failure-power"),
    "e5_pow" : qc.event("Short circuit",3.33e-4,"LOM",count = 1, source="~Estimation"),

	"e1_ls" : qc.event("Oxygen supply failure",2.48e-5,"LOC",count = 3,redundancy = 2, source="fta:launchvehicle_risk"),
    "e2_ls" : qc.event("Capsule fire",6.19e-6,"LOM",count = 1, source="fta:launchvehicle_risk"),
    "e3_ls" : qc.event("Atmospheric control failure",2.48e-5,"LOM",count = 2,redundancy = 1, source="fta:launchvehicle_risk"),
    "e4_ls" : qc.event("Capsule radiator failure",1.86e-5,"LOC",count = 2,redundancy = 1, source="fta:launchvehicle_risk"),
    "e5_ls" : qc.event("Waste management failure",4.95e-5,"LOM",count = 1, source="fta:launchvehicle_risk"),

	"e1_ab" : qc.event("Abort detection need failure",5.00e-3,"LOC",count = 2, redundancy=1, source="NASAEASA"),
    "e2_ab" : qc.event("Abort engine failure",1.09e-3,"LOM",count = 6, redundancy=1, source="fta:launchvehicle_risk"),
    "e3_ab" : qc.event("Suicide burn starting too early or late",1.10e-5,"LOM",count = 1, source="fta:utm-reliability"),
    "e4_ab" : qc.event("Too high load",1.00e-4,"LOM",count = 1, source="~Estimation"),
    "e5_ab" : qc.event("Parachute deployment failure",4.00e-3,"LOM",count = 4, redundancy=1, source="fta:parachutes-risk-mitig"),

	"e1_comp" : qc.event("Computer failure",1.67e-4,"LOV",count = 4, redundancy=3, source="midtermreport"),
    "e2_comp" : qc.event("Critical software errors",8.33e-4,"LOV",count = 1, source="midtermreport"),
    "e3_comp" : qc.event("Sudden processing errors",8.33e-4,"LOM",count = 1, source="midtermreport"),

}
compined_dict = {
    "c1_s" : qc.comp("Thrust structure failure and capsule not separating","LOC",["e2_s","e6_s"]),
	"c2_s" : qc.comp("Landing legs not deploying and capsule not separating","LOC",["e2_s","e6_s"]),
    "c3_s" : qc.comp("Pipe leakage with skirt buckling","LOV",["e3_s","e5_s"]),
	"c1_en_s" : qc.comp("Fluid leak with malfunctioning hatches","LOC",["e4_en","e13_s","e14_s"]),
	"c2_en_s" : qc.comp("Control System failure and capsule not separating","LOC",["e2_en","e6_s"]),
	"c3_en_s" : qc.comp("Engine control system failure and thrust structure failure","LOC",["e2_en","e2_s"]),
	"c4_en_s" : qc.comp("Gimbal bearing failure and thrust structure failure","LOC",["e8_en","e2_s"]),
    "c1_gnc" : qc.comp("Shutoff valve and critical failure","LOC",["e9_gnc","e1_gnc"]),
    "c2_gnc" : qc.comp("Wrong attitude/position during landing","LOV",["e2_gnc","e13_gnc"]),
    "c3_gnc" : qc.comp("Wrong attitude/position during docking","LOV",["e1_gnc","e13_gnc"]),
	"c1_en_th" : qc.comp("Propellant leak with damaged TPS","LOC",["e4_en","e1_th"]),
	"c1_gnc_th" : qc.comp("RCS propellant leak with damaged TPS","LOC",["e11_gnc","e1_th"]),
	"c1_gnc_s" : qc.comp("Guidance algorithm errors and capsule ring not deploying","LOC",["e12_gnc","e15_s"]),
	"c2_gnc_s" : qc.comp("Lading camera failure and Landing leg collapsing","LOC",["e2_gnc","e4_s"]),
	"c1_gnc_en" : qc.comp("RCS thruster failure and propellant leakage","LOC",["e7_gnc","e4_en"]),
	"c2_gnc_en" : qc.comp("Gimbal bearing failure and guidance algorithm errors","LOC",["e13_gnc","e8_en"]),
	"c1_com_ls" : qc.comp("No line of sight post abort landing and Oxygen supply failure","LOC",["e5_com","e1_ls"]),
	"c2_com_ls" : qc.comp("No line of sight post abort landing and Capsule radiator failure","LOC",["e5_com","e4_ls"]),
	"c1_com_en" : qc.comp("Short circuit and propellant leakage","LOC",["e5_pow","e4_en"]),
	"c1_ab_s" : qc.comp("Buckling of the capsule and too high loads on abort","LOC",["e4_ab","e11_s"]),
	"c1_ab_s" : qc.comp("Abort engine failure and thrust structure failure","LOC",["e2_ab","e2_s"]),
	"c1_ab_LOV" : qc.comp("Abort engine failure and LOV event","LOC",["e2_ab","all_LOV"]),
	"c1_s_LOV" : qc.comp("Capsule not separating and LOV event","LOC",["e6_s","all_LOV"]),
	

}
test = qc.PRA(event_dict,compined_dict,["LOM","LOV","LOC"])
test.gen_table("1", "1")
print((1-np.array(test.proability))*100)