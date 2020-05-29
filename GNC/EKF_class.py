import numpy as np
import sympy as smp
from numpy import linalg as npl
from sympy.solvers import solve
from sympy.utilities.lambdify import lambdify
import copy

class Jacobian:
    def __init__(self,sys):
        self.main = []
        for eq in sys.sym:
                temp = []
                for r in sys.var:
                    temp.append(lambdify(sys.var,smp.diff(eq,r)))
                self.main.append(temp)
        self.size = np.shape(self.main)
       
    def eval(self,x):
        M = np.ndarray(self.size)
        for i in range(self.size[0]):
                for j in range(self.size[1]):
                    M[i][j] = self.main[i][j](*x)
        return M

class vector_func:
    def __init__(self,function,vector):
        self.sym = function
        self.var = vector

    def subs(self,args):
        u_list = []
        for eq in self.sym:
            if type(eq) is not (int or float):
                u_list += [eq.subs(args)]
            else:
                u_list += [eq] 
        var_new = copy.deepcopy(self.var)
        for sub in args:
            var_new.remove(sub[0])
            if type(sub[1]) != (int or float):
                for el in sub[1].free_symbols:
                    if el not in var_new:
                        var_new += [el]
        return vector_func(u_list,var_new)

    def subs_func(self,func,args = 0):
        if args == 0:
            args = self.var
        u_list = []
        sub_list = [(args[i],func.sym[i]) for i in range(len(func.sym))]
        for eq in self.sym:
            if type(eq) is not (int or float):
                u_list += [eq.subs(sub_list)]
            else:
                u_list += [eq] 
        var_new = copy.deepcopy(self.var)
        for sub in sub_list:
            var_new.remove(sub[0])
            if type(sub[1]) != (int or float):
                for el in sub[1].free_symbols:
                    if el not in var_new:
                        var_new += [el]
        return vector_func(u_list,var_new)

    def lamb(self,ref_list = 0):
        if ref_list == 0:
            ref_list = self.var
        self.main = []
        for eq in self.sym:
            self.main.append(lambdify(ref_list,eq))
        return 0
    
    def eval(self,x):
        return np.array([eq(*x) for eq in self.main])
            
class EKF:
    def __init__(self,f,h,Q,R,x0):
        self.x_var = f.var
        self.x = x0
        self.f = f
        self.f.lamb()
        self.h = h
        self.h.lamb()
        self.P = np.ndarray([len(f.main),len(f.main)])
        self.Q = Q
        self.J = Jacobian(f)
        self.R = R
        self.n = len(f.main)

    def step(self):
        Jn  = self.J.eval(self.x)
        self.x = self.f.eval(self.x)
        self.P = Jn*self.P*np.transpose(Jn) + self.Q

    def corr(self,z):
        Jn = self.J.eval(self.x)
        #print("Jn")
        #print(Jn)

        self.K = self.P*Jn*npl.inv(Jn*self.P*np.transpose(Jn)+self.R)
        #print("K")
        #print(self.K)
        self.Pn = (np.eye(self.n)-self.K*Jn)*self.P
        #print("Pn")
        #print(self.Pn)
        self.P = self.Pn
        #print("x")
        #print(self.x)
        self.x = self.x + np.dot(self.K,(z-self.h.eval(self.x)))