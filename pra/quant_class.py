import importlib
import os.pat
import itertools
import shutil
class event:
    def __init__(self,description,proability,consqequence,count,redundancy = 0)
        self.desc = description
        self.prob = proability
        self.con = consqequence
        self.count = count
        self.con_n = 0
        self.red = redundancy
    class comp:
    def __init__(self,description,consqequence,event_list)
        self.desc = description
        self.con = consqequence
        self.event_list = event_list

class PRA:
    def __init__(self,evnt_dict,com_dict,con_namelist):
        self.con_list = con_namelist
        self.con_n = range(len(namelist))
        self.proability = [0]*len(namelist)
        self.e_dict = evnt_dict
        self.c_dict = com_dict
        for key in self.e_dict:
            event = self.e_dict[key]
            for n in self.con_n:
                if event.con == self.con_list[n]:
                    event.con_n = n
                    break:
        for i in self.con_n:
            pr = 1
            for key in self.e_dict:
                event = self.e_dict[key]
                if event.con_n == i:
                    temp = 1
                    for n in range(1+event.redundancy):
                        pr *= (1-(1-event.prob)**(event.count-n))
            for key in c_dict = 



        



        

                    




    
