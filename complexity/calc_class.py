import numpy as np

i=0
class line:
    def __init__(self,complexity,direction="uni"):
        self.cm = complexity
        if direction=="uni":
            self.lin = 1
        elif direction == "bi":
            self.lin = 2
        else:
            raise Exception("not valid input")

class system:
    def __init__(self,complexity):
        global i
        self.connections_complexity = []
        self.connections_destination = []
        self.complexity = complexity 
        self.id = i
        i += 1
        

    def add_c(self,conection,destanation):
        self.connections_complexity.append(conection.cm*conection.lin)
        destanation.connections_complexity.append(conection.cm*conection.lin)
        self.connections_destination.append(destanation.id)
        destanation.connections_destination.append(self.id)


class complexity:
    def __init__(self,l):
        M = np.zeros((len(l),len(l)))
        for i in l:
            for k in i.connections_destination:
                M[i.id,k]=1
        vec = np.linalg.eigvals(M)
        self.E = 0
        for val in vec:
            self.E += abs(val)
        self.C1 = 0
        for i in l:
            self.C1 += i.complexity
        self.C2 = 0
        for i in l:
            for val in i.connections_complexity:
                self.C2 += val/2
        self.structural = self.C1 +self.C2*self.E/len(l)