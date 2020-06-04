import pandas as pd
import numpy as np

class dataset:
    def __init__(self):
        acceldata = pd.read_csv("../astrodynamics/accelerations.csv")
        self.t = np.array(acceldata["time"])
        self.r_p = np.transpose(np.array([acceldata["x"],acceldata["y"],acceldata["z"]]))
        lat = -41*np.pi/180
        lon = 23.5*np.pi/180
        R0 = np.array([3380.9,0,0])
        ry = np.array([[np.cos(lat),0,np.sin(lat)],
                       [0,1,0],
                       [-np.sin(lat),0,np.cos(lat)]])
        rz = np.array([[np.cos(lon),-np.sin(lon),0],
                       [np.sin(lon),np.cos(lon),0],
                       [0,0,1]])    
        mars_r = np.dot(np.dot(rz,ry),R0)
        acc = np.array([acceldata["downward a"],acceldata["forwards a"]])
        self.acc = np.array([-acc[0][i]*self.r_p[i]/(np.linalg.norm(self.r_p[i])*10**3) + (mars_r-self.r_p[i])*acc[1][i]//(np.linalg.norm(mars_r-self.r_p[i])*10**3) for i in range(len(self.t))])