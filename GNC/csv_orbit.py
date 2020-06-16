import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class dataset:
    def __init__(self):
        acceldata = pd.read_csv("../astrodynamics/from_orbit.csv")
        self.t = np.array(acceldata["time"])
        self.dt = self.t[1]-self.t[0]
        self.r_p = np.transpose(np.array([acceldata["x"],acceldata["y"],acceldata["z"]]))/10**3
        acc_test = np.array([(self.r_p[i+1]-2*self.r_p[i]+self.r_p[i-1])/self.dt**2 for i in range(1,len(self.r_p)-1)])
        a0 = [(2*self.r_p[0]-5*self.r_p[1]+4*self.r_p[2]-1*self.r_p[3])/self.dt**2]
        af = [(2*self.r_p[-1]-5*self.r_p[-2]+4*self.r_p[-3]-1*self.r_p[-4])/self.dt**2]
        self.acc = np.transpose(np.concatenate([a0,acc_test,af]))
        lat = -42.5*np.pi/180
        lon = 25.5*np.pi/180
        self.R0 = np.array([3380.9,0,0])
        ry = np.array([[np.cos(lat),0,np.sin(lat)],
                       [0,1,0],
                       [-np.sin(lat),0,np.cos(lat)]])
        rz = np.array([[np.cos(lon),-np.sin(lon),0],
                       [np.sin(lon),np.cos(lon),0],
                       [0,0,1]])
        self.trans_mat= np.dot(rz,ry)
        self.mars_r = np.dot(self.trans_mat,self.R0)
        self.v0 = (-3*self.r_p[0]+4*self.r_p[1]-self.r_p[2])/(2*self.dt)
        self.a0 = np.transpose(self.acc)[0]