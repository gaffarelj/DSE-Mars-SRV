from decimal import *
import numpy as np
import copy

class event:
    def __init__(self,description,proability,consqequence,count=1,redundancy = 0,source=""):
        self.desc = description
        self.count = count
        self.red = redundancy
        self.prob_single = proability
        pr = 1
        for i in range(redundancy+1):
            pr *= (1-(1-self.prob_single)**(self.count-i))
        self.prob = 1 - pr
        self.con = consqequence
        self.con_n = 0
        self.source = source
    def redo_prob(self):
        pr = 1
        for i in range(self.red+1):
            pr *= (1-(1-self.prob_single)**(self.count-i))
        self.prob = 1 - pr

class comp:
    def __init__(self,description,consqequence,event_list):
        self.desc = description
        self.con = consqequence
        self.event_list = event_list
        self.con_n = 0

        self.val_comp = False
        self.prob_val = 0
    
    def prob(self, pra):
        if not(self.val_comp):
            pr = 1
            for key in self.event_list:
                if key in pra.e_dict:
                    pr *= 1- pra.e_dict[key].prob
                elif key in pra.c_dict:
                    if pra.c_dict[key].val_comp:
                        pr *= 1-pra.c_dict[key].prob_val     
                    else:
                        pr *= 1-pra.c_dict[key].prob(pra)
                elif key[:4] == "all_":
                    if key[4:] in pra.con_list:
                        pr *= 1-pra.proability[pra.con_list.index(key[4:])]
            self.prob_val = 1- pr
            self.val_comp = True
            #print(self.prob_val)
        return self.prob_val
    
    def name_func(self, pra):
        self.name = "("
        for key in self.event_list:
            if key in pra.e_dict:
                self.name += key.replace("_", "$_{\\text{") + "}}$"
            elif key in pra.c_dict:
                self.name += pra.c_dict[key].name_func(pra)
            elif key[:4] == "all_":
                self.name += key.replace("_", " ")
            if key != self.event_list[-1]:
                self.name += " AND "
        self.name += ")"
        return self.name



class PRA:
    def __init__(self,evnt_dict,com_dict,con_namelist):
        self.con_list = con_namelist
        self.con_n = range(len(con_namelist))
        self.proability = [0]*len(con_namelist)
        self.e_dict = evnt_dict
        self.c_dict = com_dict
        for key in self.e_dict:
            event = self.e_dict[key]
            for n in self.con_n:
                if event.con == self.con_list[n]:
                    event.con_n = n
                    break

        for key in self.c_dict:
            comp = self.c_dict[key]
            for n in self.con_n:
                if comp.con == self.con_list[n]:
                    comp.con_n = n
                    break

        for i in self.con_n:
            pr = 1
            for key in self.e_dict:
                event = self.e_dict[key]
                if event.con_n == i:
                    pr *= event.prob
            for key in self.c_dict:
                comp = self.c_dict[key]
                
                if comp.con_n == i:
                    #print(key)
                    pr *= comp.prob(self)
            self.proability[i] = pr

    def sens(self,n, div):
        temp_ref_list = []
        for i in range(n):
            temp_e_dict = copy.deepcopy(self.e_dict)
            for key in temp_e_dict:
                temp_e_dict[key].prob_single = self.e_dict[key].prob_single + np.random.normal(0,div*self.e_dict[key].prob_single)
                temp_e_dict[key].redo_prob()
            temp_pra = PRA(temp_e_dict,self.c_dict,self.con_list)
            temp_ref_list.append(temp_pra.proability)
        self.ref_list = np.array(temp_ref_list)
        self.std = np.std(self.ref_list,axis = 0)
        self.mu = np.average(self.ref_list,axis = 0)


    def gen_table(self,caption, label):
        print("\\begin{longtable}[H]{|l|l|l|l|}")
        print("\\caption{" + caption + "}")
        print("\\label{tab:"+label+"} \\\\")
        print("\\hline")
        for i in range(len(self.con_list)):
            con = self.con_list[i]
            print("\\multicolumn{4}{|c|}{\\textbf{" + con + "}}\\\\\\hline")
            print("\\textbf{ID} & \\textbf{Source} & \\textbf{Description} &  \\textbf{Probability [-]} \\\\\\hline")
            for key in self.e_dict:
                event = self.e_dict[key]
                if event.con == con:
                    if event.source != "":
                        if event.source[0] == "~":
                            src = event.source[1:]
                        else:
                            src = "\\cite{"+ event.source +"}"
                    else:
                        src = ""
                    if event.prob == 1:
                        print(key.replace("_", "$_{\\text{") + "}}$ & " + src + " & " + event.desc + " & " + "0E0" + "\\\\\\hline")
                    else:
                        print(key.replace("_", "$_{\\text{") + "}}$ & " + src + " & " + event.desc + " & " + "{:.2E}".format(Decimal((1-event.prob))) + "\\\\\\hline")
            
            for key in self.c_dict:
                comp = self.c_dict[key]
                if comp.con == con:
                    if comp.prob_val == 1:
                        print("\\multicolumn{2}{|l|}{" + comp.name_func(self) + "} & " + comp.desc + " & " + "0E0" + "\\\\\\hline")
                    else:
                        print("\\multicolumn{2}{|l|}{" + comp.name_func(self) + "} & " + comp.desc + " & " + "{:.2E}".format(Decimal((1-comp.prob_val))) + "\\\\\\hline")
            print("\\multicolumn{3}{|l|}{\\textbf{Total " + con + " }} & \\textbf{" + "{:.2E}".format(Decimal((1-self.proability[i]))) + "} \\\\\\hline")
        print("\end{longtable}")





        

                    




    
