
class event:
    def __init__(self,description,proability,consqequence,count,redundancy = 0):
        self.desc = description
        self.count = count
        pr = 1
        for i in range(redundancy+1):
            pr *= (1-(1-proability)**(self.count-i))
        self.prob = 1 - pr
        self.con = consqequence
        
        self.con_n = 0
        self.red = redundancy

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
            print(self.prob_val)
        return self.prob_val
    
    def name_func(self, pra):
        self.name = "("
        for key in self.event_list:
            if key in pra.e_dict:
                self.name += key
            elif key in pra.c_dict:
                self.name += pra.c_dict[key].name_func(pra)
            elif key[:4] == "all_":
                self.name += key
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
                    print(key)
                    pr *= comp.prob(self)
            self.proability[i] = pr
    def gen_table(self,caption, label):
        print("\\begin{table}[H]")
        print("\\centering")
        print("\\caption{" + caption + "}")
        print("\\begin{tabular}{|l|c|c|} \\hline")
        for i in range(len(self.con_list)):
            con = self.con_list[i]
            print("\\multicolumn{2}{|c|}{\\textbf{" + con + "}} & " + str(round(1/(1-self.proability[i]))) + "\\\\\\hline")
            for key in self.e_dict:
                event = self.e_dict[key]
                if event.con == con:
                    print(key + " & " + event.desc + " & " + str(round(1/(1-event.prob))) + "\\\\\\hline")
            
            for key in self.c_dict:
                comp = self.c_dict[key]
                if comp.con == con:
                    print(comp.name_func(self) + " & " + comp.desc + " & " + str(round(1/(1-comp.prob_val))) + "\\\\\\hline")
        print("\\label{tab:"+label+"}")
        print("\end{tabular}")
        print("\end{table}")





        

                    




    
