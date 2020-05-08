import numpy as np


class color:
	def __init__(self, code, name):
		self.HTML = code
		self.name = name
	
class param:
	def __init__(self, name, weight, func="LRTS",  direc="HB",  p=1):
		self.name = name
		self.func = func
		self.dir = direc
		self.p = p
		self.weight = weight

	def func_eval(self, Hv, Lv, evalv):
		if self.dir == "LB":
			Hv, Lv = Lv, Hv
		if self.func == "LRTS":
			return (evalv - Lv) / (Hv - Lv)
		elif self.func == "IRTS":
			return (1 - exp(-(evalv - Lv) / self.p)) / (1 - exp(-(Hv - Lv) / self.p))
		elif self.func == "DRTS":
			return (1 - exp(-(Hv - evalv) / self.p)) / (1 - exp(-(Hv - Lv) / self.p))
		else:
			raise Exception("not valid input")
	
	def set_colors(self, eval_list, color_list):
		self.values = eval_list
		self.color = []
		for val in eval_list:
			if val != 1:
				self.color.append(color_list[int(val / len(color_list))])
			else:
				self.color.append(color_list[int(val / len(color_list)) - 1])
		

class design:
	def __init__(self, name, sourcelist):
		self.name = name
		self.sourcelist = sourcelist


class tradeoff:
	def __init__(self, design_list, param_list, out="python", width=8, color_list=[]):
		for i in range(len(param_list)):
			param = param_list[i]
			param_val = []
			for design in design_list:
				param_val.append(design.sourcelist[i])
			Lv = min(param_val)
			Hv = max(param_val)
			eval_list = []
			for val in param_val:
				eval_list.append(round(param.func_eval(Hv, Lv, val),  5))

			param.set_colors(eval_list, color_list)
			#param.values = eval_list

			if out == "python":
				output = param.name + ",  Actual value:  "
				for val in param_val:
					output += str(val) + ",    "
				print(output)
				output = param.name + ",  scaled and weighted value:  "
				for val in eval_list:
					output += str(val) + ",    "
				print(output)
			
		if out == "latex":
			print("\\begin{table}[]")
			print("\caption{}")
			print("\label{tab:my-table}")
			print("\\begin{adjustbox}{width=\linewidth, center}")

			output = "\\begin{tabular}{|c|l|"
			for param in param_list:
				output += "p{" + str(width * param.weight) + "cm}|"
				output += "p{" + str(width * param.weight) + "cm}|"
			output +="c|}\hline"
			print(output)

			output = "\multicolumn{2}{|c|}{\\textbf{Criteria}}"
			for param in param_list:
				output += "& \multicolumn{2}{c|}{}"
			print(str(output) + "\\\\")
			output = "\cline{1-2}\multicolumn{2}{|l|}{\\textbf{Design Option}}"
			for param in param_list:
				output += "& \multicolumn{2}{c|}{\multirow{-2}{*}{" + param.name + "}}"
			print(str(output) + "& \multirow{-2}{*}{\\textbf{Total}} \\\\ \hline")
			for i in range(len(design_list)):
				design = design_list[i]
				output = "\multicolumn{2}{|c|}{}"
				end_output = ""
				k = 4
				for param in param_list:
					output += "   & \cellcolor[HTML]{" + str(param.color[i].HTML) + "} & \cellcolor[HTML]{" + str(param.color[i].HTML) + "}" + str(param.color[i].name) + ""
					end_output += " \cline{" + str(k) + "-" + str(k) + "} "
					k += 2
				print(str(output) + " & \\\\" + str(end_output))
				output = "\multicolumn{2}{|c|}{}"
				for param in param_list:
					output += "   & \multicolumn{2}{c|}{\cellcolor[HTML]{" + str(param.color[i].HTML) + "}}"
				print(str(output) + "& \\\\")
				output = "\multicolumn{2}{|c|}{\multirow{-3}{*}{" + str(design.name) + "}}"
				total = 0
				for param in param_list:
					output += "   &\multicolumn{2}{c|}{\multirow{-2}{*}{\cellcolor[HTML]{" + str(param.color[i].HTML) + "}" + str(param.values[i])[:5] + "}}"
					total += param.values[i] * param.weight
				print(str(output) + " & \multirow{-3}{*}{" + str(total)[:5] + "} \\\\ \hline")
			print("\end{tabular}")
			print("\end{adjustbox}")
			print("\end{table}")