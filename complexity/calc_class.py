"Written By Lasse Landergren"
import numpy as np

l = []
def reset():
	global l
	l = []

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
		self.connections_complexity = []
		self.connections_destination = []
		self.complexity = complexity
		self.id = len(l)
		l.append(self)
		
	def add_c(self,conection,destinations):
		if type(destinations) == system:
			destinations = [destinations]
		if type(destinations) != list:
			raise Exception("destination invalid")
		for destination in destinations:
			self.connections_complexity.append(conection.cm*conection.lin)
			destination.connections_complexity.append(conection.cm*conection.lin)
			self.connections_destination.append(destination.id)
			destination.connections_destination.append(self.id)


class complexity:
	def __init__(self):
		M = np.zeros((len(l),len(l)))
		for i in l:
			for k in i.connections_destination:
				M[i.id, k] = 1
		vec = np.linalg.eigvals(M)
		self.E = sum(abs(vec))
		self.C1 = sum([i.complexity for i in l])
		self.C2 = 0
		for i in l:
			for val in i.connections_complexity:
				self.C2 += val*i.complexity/2
		self.structural = self.C1 + self.C2*self.E/len(l)
		self.average = self.structural/len(l)