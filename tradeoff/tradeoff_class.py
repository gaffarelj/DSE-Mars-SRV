import numpy as np
import copy
import multiprocessing as mp

class color:
	def __init__(self, code, name):
		self.HTML = code
		self.name = name
	
class param:
	def __init__(self, name, weight, func="LRTS",  direc="HB",  p=1,Limitype ="minmax"):
		self.name = name
		self.func = func
		self.dir = direc
		self.p = p
		self.weight = weight
		self.Ltype = Limitype
		self.val_in = []
		self.val_out = []
	
	def stat(self):
		self.sd = np.std(self.val_in)
		self.mu = np.average(self.val_in)

		if self.Ltype == "minmax":
			self.Lv, self.Hv = min(self.val_in), max(self.val_in)
		elif self.Ltype == "1SD":
			self.Lv, self.Hv = self.mu-self.sd, self.mu+self.sd
		elif self.Ltype == "2SD":
			self.Lv, self.Hv = self.mu-2*self.sd, self.mu+2*self.sd
		elif self.Ltype == "3SD":
			self.Lv, self.Hv = self.mu-3*self.sd, self.mu+3*self.sd
		else:
			raise Exception("not valid boundary determination method")



	def func_eval(self, evalv):
		if (evalv <= self.Lv and self.dir == "HB") or (evalv >= self.Hv and self.dir == "LB") :
			return 0
		if (evalv >= self.Hv and self.dir == "HB") or (evalv <= self.Lv and self.dir == "LB") :
			return 1

		if self.dir == "LB":
			temHv, temLv = self.Lv, self.Hv
		else:
			temHv, temLv = self.Hv, self.Lv

		if self.func == "LRTS":
			return (evalv - temLv) / (temHv - temLv)
		elif self.func == "IRTS":
			return (1 - np.exp(-(evalv - temLv) / self.p)) / (1 - np.exp(-(temHv - temLv) / self.p))
		elif self.func == "DRTS":

			return (1 - np.exp(-(temHv - evalv) / self.p)) / (1 - np.exp(-(temHv- temLv) / self.p))
		else:
			raise Exception("not valid scoring scheme")
	
	def set_colors(self, color_list):
		self.color = []
		for val in self.val_out:
			if val != 1:
				self.color.append(color_list[int(val * len(color_list))])
			else:
				self.color.append(color_list[int(val * len(color_list)) - 1])
	
	
class design:
	def __init__(self, name, sourcelist):
		self.name = name
		self.sourcelist = sourcelist


class tradeoff:
	def __init__(self, design_list, param_list):
		self.param_list = param_list
		self.design_list = design_list


	def get_tradeoff(self):
		self.total = np.zeros(len(self.design_list))
		for i in range(len(self.param_list)):
			param = self.param_list[i]
			param.val_in = np.array([design.sourcelist[i] for design in self.design_list])
			param.stat()
			param.val_out = np.array([param.func_eval(val) for val in param.val_in])
			self.total += param.val_out*param.weight
		
	
	def get_output(self,language = "python",color_list=[],width=8):
		if language == "python":
			for param in self.param_list:
				print(param.name, ",\t actual value:", end="\t", sep="")
				for val in param.val_in:
					print(val, end=",\t")
				print()
				print(param.name, ",\t scaled value:", end="\t", sep="")
				for val in param.val_out:
					print(round(val,5), end=",\t")
				print()
			print("\t final value:", end="\t", sep="")
			for val in self.total:
				print(round(val,5), end=",\t")
			print()
		if language == "latex":
			if len(color_list)==0:
				raise Exception("color_list is mandatory for Latex output")
			for param in self.param_list:
				param.set_colors(color_list)
			print("\\begin{table}[]")
			print("\caption{}")
			print("\label{tab:my-table}")
			print("\\begin{adjustbox}{width=\linewidth, center}")

			output = "\\begin{tabular}{|c|l|"
			for param in self.param_list:
				output += "p{" + str(width * param.weight) + "cm}|"
				output += "p{" + str(width * param.weight) + "cm}|"
			output +="c|}\hline"
			print(output)

			output = "\multicolumn{2}{|c|}{\\textbf{Criteria}}"
			for param in self.param_list:
				output += "& \multicolumn{2}{c|}{"+ param.name +"}"
			print(str(output) + "&\\\\ \cline{1-2}")
			output = "\multicolumn{2}{|l|}{\\textbf{Design Option}}"
			for param in self.param_list:
				output += "& \multicolumn{2}{c|}{\\textit{("+ str(round(param.Lv,2)) + " / " + str(round(param.Hv,2)) + ")," 
				if param.dir == "HB":
					output += " High Best}}"
				else:
					output += " High Best}}"
				
			print(str(output) + "& \multirow{-2}{*}{\\textbf{Total}} \\\\ \hline")
			for i in range(len(self.design_list)):
				design = self.design_list[i]
				output = "\multicolumn{2}{|c|}{}"
				end_output = ""
				k = 4
				for param in self.param_list:
					output += "   & \cellcolor[HTML]{" + str(param.color[i].HTML) + "} & \cellcolor[HTML]{" + str(param.color[i].HTML) + "}" + str(param.color[i].name) + ""
					end_output += " \cline{" + str(k) + "-" + str(k) + "} "
					k += 2
				print(str(output) + " & \\\\" + str(end_output))
				output = "\multicolumn{2}{|c|}{}"
				for param in self.param_list:
					output += "   & \multicolumn{2}{c|}{\cellcolor[HTML]{" + str(param.color[i].HTML) + "}}"
				print(str(output) + "& \\\\")
				output = "\multicolumn{2}{|c|}{\multirow{-3}{*}{" + str(design.name) + "}}"
				for param in self.param_list:
					output += "   &\multicolumn{2}{c|}{\multirow{-2}{*}{\cellcolor[HTML]{" + str(param.color[i].HTML) + "}" + str(round(param.val_in[i])) + " / " + str(round(param.val_out[i],5)) + "}}"
				print(str(output) + " & \multirow{-3}{*}{" + str(round(self.total[i],5)) + "} \\\\ \hline")
			print("\end{tabular}")
			print("\end{adjustbox}")
			print("\end{table}")

class sensitivity:
	def __init__(self,tradeoff,samples = 10000):
		self.tro = tradeoff
		self.n = samples
		self.to_tech = False
		self.to_p = False
		self.to_weights = False

	def addto_technical(self,variation):
		self.to_tech = True
		self.to_tech_var = variation

	def addto_p(self,variation):
		self.to_p = True
		self.to_p_var = variation

	def addto_weights(self,variation):
		self.to_weights = True
		weight_list = np.array([param.weight for param in self.tro.param_list])
		self.to_weights_var = variation*np.std(weight_list)

	def sens(self,n):
			tro_temp = copy.deepcopy(self.tro)
			if self.to_p:
				for param in tro_temp.param_list:
					param.p = np.random.normal(param.p,self,self.to_p_var)

			if self.to_weights:
				total = 0
				for param in tro_temp.param_list:
					param.weight = np.random.normal(param.weight,self.to_weights_var)
					total += param.weight
				
				for param in tro_temp.param_list:
					param.weight /= total

			if self.to_tech:
				for design in tro_temp.design_list:
					for i in range(len(design.sourcelist)):
						design.sourcelist[i] = np.random.normal(design.sourcelist[i],self.tro.param_list[i].sd*self.to_tech_var)
			
			tro_temp.get_tradeoff()
			ret = np.zeros(len(tro_temp.design_list))
			ret[np.where(tro_temp.total == np.amax(tro_temp.total))] = 1
			return ret

	def get_sens(self):
		pool = mp.Pool(mp.cpu_count())
		self.per = np.sum(pool.map(self.sens,range(self.n)),0)
		self.per /= self.n
	
	def get_RMS(self):
		self.RMS = np.zeros(len(self.tro.design_list))
		for param in self.tro.param_list:
			self.RMS += np.multiply(param.val_in-param.mu,param.val_in-param.mu)/(param.sd*param.sd)
		

