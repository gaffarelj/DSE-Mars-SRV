
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
    
    def prob(self, c_dict, e_dict):
        if not(self.val_comp):
            pr = 1
            for key in self.event_list:
                if key in e_dict:
                    pr *= 1- e_dict[key].prob
                elif key in c_dict:
                    if c_dict[key].val_comp:
                        pr *= 1-prob_val
                        
                    else:
                        pr *= 1-c_dict[key].prob(c_dict,e_dict)
            self.prob_val = 1- pr
            self.val_comp = True
            print(self.prob_val)
        return self.prob_val
                



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
                    pr *= comp.prob(self.c_dict,self.e_dict)
            self.proability[i] = pr


        



        

                    




    
