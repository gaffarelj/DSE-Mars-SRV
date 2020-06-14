# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:02:44 2020

@author: EnterGamer
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

g = 9.81
ax= 9

F = 32000000
CG =8.04
##---Mass elements---###
R= 3
A1_Engines = 3150
A2_Engines_2 = 720
A3_Thrust_Str = 515
A4_Skirt = 1000
A5_RCS = 400
A6_Tank_f = 1100  + 200 +450
A7_Tank_o = 2500 + 601 +300
A8_Skirt_2 = 1000 + 298.8 
A9_Adapter = 600
A10_Hydra = 1000
A11_TPS_ch = 1400
A12_Landing_legs = 3600
A13_Capsule1 = 4768.11
A14_Capsule2 = 7110.688
A15_Capsule3 = 2393.5794
B13_Fuel = 20171.83164
B14_Oxi = 76652.96022

M = A1_Engines+A2_Engines_2+A3_Thrust_Str+A4_Skirt+A5_RCS+A6_Tank_f+A7_Tank_o+A8_Skirt_2+A9_Adapter+A10_Hydra+A11_TPS_ch+A12_Landing_legs+A13_Capsule1+A14_Capsule2+A15_Capsule3
print(M)
M_f = M+B13_Fuel+B14_Oxi

A1_s = 0
A1_f = 3

A2_s = 2.5
A2_f = 3

A3_s = 2.5
A3_f = 3

A4_s = 2.5
A4_f = 5.5

A5_s = 3
A5_f = 3.3

A6_s = 3.3
A6_f = 5.5

A7_s = 5.5
A7_f = 11.3

A8_s = 9.1
A8_f = 12.1

A9_s = 11.8
A9_f = 12.3

A10_s = 11.3
A10_f = 12.1

A11_s = 2.5
A11_f = 12.3

A12_s = 2.5
A12_f = 8.5

A13_s = 12.3
A13_f = 13.5

A14_s = 13.5
A14_f = 15.2

A15_s = 15.2
A15_f = 17.8225

B13_s = 3.3
B13_f = 5.5

B14_s = 5.5
B14_f = 11.3

p_1 = A1_Engines/ (A1_f - A1_s)
p_2 = A2_Engines_2 / (A2_f - A2_s)
p_3 = A3_Thrust_Str / (A3_f - A3_s)
p_4 = A4_Skirt / (A4_f - A4_s)
p_5 = A5_RCS / (A5_f - A5_s)
p_6 = A6_Tank_f / (A6_f - A6_s)
p_7 = A7_Tank_o / (A7_f - A7_s)
p_8 = A8_Skirt_2 / (A8_f - A8_s)
p_9 = A9_Adapter / (A9_f - A9_s)
p_10 = A10_Hydra / (A10_f - A10_s)
p_11 = A11_TPS_ch / (A11_f - A11_s)
p_12 = A12_Landing_legs / (A12_f - A12_s)
p_13 = (B13_Fuel / (B13_f - B13_s))
p_14 = (B14_Oxi / (B14_f - B14_s))
p_15 = (A13_Capsule1/ (A13_f-A13_s))
p_16 = (A14_Capsule2/ (A14_f - A14_s))
p_17 = (A15_Capsule3/ (A15_f-A15_s))


pz_1 = p_1
pz_2 = (p_1+p_2+p_3+p_4+p_11+p_12)
pz_3 = (p_4+p_5+p_11+p_12)
pz_4 = (p_4+p_6+p_11+p_12)
pz_5 = (p_7 + p_11 + p_12)
pz_6 = (p_7 + p_11)
pz_7 = (p_7 + p_8 + p_11)
pz_8 = (p_8 + p_10 + p_11)
pz_9 = (p_8 + p_9 + p_10 + p_11)
pz_10 = (p_9 + p_11)

z_1 = 2.5
z_2 = 0.5
z_3 = 0.3
z_4 = 2.2
z_5 = 3.0
z_6 = 0.6
z_7 = 2.2
z_8 = 0.5
z_9 = 0.3
z_10 = 0.2
z_11 = 5.8
z_12 = 1.2
z_13 = 1.7
z_14 = 2.6225

dz = 0.01

z_l = np.arange(0, 17.83, dz)

p_zf = []
for z in z_l:
    if z >= 0 and z< 2.5:
        p_zf.append(pz_1)
    if z >= 2.5 and z< 3:
        p_zf.append(pz_2)
    if z >= 3 and z< 3.3:
        p_zf.append(pz_3)
    if z >= 3.3 and z< 5.5:
        p_zf.append(pz_4+p_13)
    if z >=5.5 and z < 8.5:
        p_zf.append(pz_5+p_14)
    if z >=8.5 and z < 9.1:
        p_zf.append(pz_6+p_14)
    if z >= 9.1 and z< 11.3:
        p_zf.append(pz_7+p_14)
    if z >= 11.3 and z< 11.8:
        p_zf.append(pz_8)
    if z >= 11.8 and z < 12.1:
        p_zf.append(pz_9)
    if z>= 12.1 and z < 12.3:
        p_zf.append(pz_10)
    if z>= 12.3 and z < 13.5:
        p_zf.append(p_15)
    if z>= 13.5 and z < 15.2:
        p_zf.append(p_16)
    if z>= 15.2 and z < 17.8225:
        p_zf.append(p_17)

alpha = 3
print(np.shape(z_1))   
V = [0]*len(p_zf)

for i in range(len(V)):
    V[i] = alpha*p_zf[i]*(z_l[i]-CG)*dz+V[i-1]




fig, ax = plt.subplots()

ax.plot(z_l, [alpha*p_zf[i]*(z_l[i]-CG) for i in range(len(z_l))])

ax.set(xlabel='z [m]', ylabel='Force [kg/m]',
       title='Moment distribution in Charon')
ax.grid()
"""
V = []
Z = []
for i in range(1,(len(tau)-1)):
    dm_dz = (tau[i+1]-tau[i-1])/(2*0.01)
    V.append(dm_dz)
    Z.append(z[i])
    

fig, ax = plt.subplots()

ax.plot(Z, V)

ax.set(xlabel='z [m]', ylabel='Force [kg/m]',
       title='Shear distribution in Charon')
ax.grid()
"""
plt.show()

