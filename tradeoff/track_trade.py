import tradeoff_class as tc

avi = tc.param(name="Availability",weight=1/3,Limitype ="fixed",Limit_val=[1,5])
acc = tc.param(name="Accuracy",weight=1/3,Limitype ="fixed",Limit_val=[1,5])
sam = tc.param(name="sample rate",weight=1/3,Limitype ="fixed",Limit_val=[1,5])
gyro = tc.design(name="Gyroscope",sourcelist=[5,3,5])
sun = tc.design(name="Sunsensor",sourcelist=[2,5,3])
star = tc.design(name="Startracker",sourcelist=[3,5,3])
maf = tc.design(name="Magnometer",sourcelist=[5,2,5])

opt = tc.design(name="Optical",sourcelist=[2,5,3])
accm = tc.design(name="Accelerometer",sourcelist=[5,3,5])

radar = tc.design(name="Radar",sourcelist=[5,5,2])
inter = tc.design(name="Radio Interforemter",sourcelist=[5,3,5])
gropt = tc.design(name="Optical (Ground)",sourcelist=[2,4,3])
lidar = tc.design(name="Lidar",sourcelist=[3,5,2])
doppler = tc.design(name="Doppler Radar",sourcelist=[5,2,5])


colors = [tc.color("EF5350", "red"), tc.color("FB8C00", "orange"), tc.color("FFEB3B", "yellow"), tc.color("8BC34A", "green"), tc.color("00BCD4", "blue")]

tradeoff_att =tc.tradeoff(design_list = [gyro,sun,star,maf],param_list= [avi,acc,sam])

tradeoff_att.get_tradeoff()
tradeoff_att.get_output(language="latex",color_list=colors,width=10,rot="hor",caption="Tradeoff Rotational sensors")
#tradeoff_att.get_output()

tradeoff_in =tc.tradeoff(design_list = [opt,accm],param_list= [avi,acc,sam])

tradeoff_in.get_tradeoff()
tradeoff_in.get_output(language="latex",color_list=colors,width=10,rot="hor",caption="Tradeoff Inertial sensors")
#tradeoff_in.get_output()

tradeoff_gr =tc.tradeoff(design_list = [radar,inter,gropt,lidar,doppler],param_list= [avi,acc,sam])

tradeoff_gr.get_tradeoff()
tradeoff_gr.get_output(language="latex",color_list=colors,width=10,rot="hor",caption="Tradeoff Ground based sensors")
#tradeoff_gr.get_output()




