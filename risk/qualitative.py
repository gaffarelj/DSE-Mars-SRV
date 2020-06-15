import csv
import os


def clean_tuple(*in_tup):
	check = [1 if (t == "" or t.isspace()) else 0 for t in in_tup]
	if 1 in check:
		return tuple([None]*len(in_tup))
	return in_tup

class risk():
	def __init__(self, id, event, conseq, prob, impact, prob_mit=(None, None), impact_mit=(None, None)):
		self.id = id
		self.event = event.strip()
		self.conseq = conseq.strip()
		self.prob = int(prob)
		self.impact = int(impact)
		self.prob_mit_exp, self.prob_mit = prob_mit
		self.impact_mit_exp, self.impact_mit = impact_mit
		self.prob_mit = self.prob if self.prob_mit is None else int(self.prob_mit)
		self.impact_mit = self.impact if self.impact_mit is None else int(self.impact_mit)

	def __str__(self):
		s =  f"Risk {self.id}: {self.event}\n"
		s += f" - Proba.: {self.prob} -> {self.prob_mit}\n"
		s += f" - Impact: {self.impact} -> {self.impact_mit}"
		return s

	def __repr__(self):
		return f"risk {self.id} ({self.prob}->{self.prob_mit}/{self.impact}->{self.impact_mit})"

class risks_list():
	def __init__(self, risk_list=[], ss_codename="xxx", ss_name="xxx", use_mitig=False, latex_file=True):
		self.risk_list = risk_list
		self.path = os.path.dirname(os.path.realpath(__file__))
		self.ss_codename = ss_codename.upper()
		self.ss_name = ss_name.title()
		self.use_mitig = use_mitig
		self.latex_file = latex_file
		self.latex = []

	def get_by_impact(self, impact):
		t = []
		for r in self.risk_list:
			if self.use_mitig:
				if r.impact_mit == impact:
					t.append(r)
			else:
				if r.impact == impact:
					t.append(r)
		return t

	def get_risks(self, fname="Example-COPY-ONLY"):
		read = csv.reader(open(f"{self.path}\input\{fname}.csv"))
		for i, row in enumerate(read):
			if i > 0:
				prob_mitig, impact_mitig = clean_tuple(row[4], row[6]), clean_tuple(row[5], row[7])
				if len(row[0]) > 0:
					r = risk(i, row[0], row[1], row[2], row[3], prob_mitig, impact_mitig)
					self.risk_list.append(r)
		self.risks_def()
		self.map()
		self.mitig(True)
		self.map()	
		self.combine_latex()

	def risks_def(self):
		ret = "\\begin{itemize}\n"
		for r in self.risk_list:
			ret += f"\t \\item \\textbf{{SRV-RISK-{self.ss_codename}-{r.id}}} {r.event.capitalize()}, results in {r.conseq.lower()}.\n"
			if r.impact_mit != r.impact or r.prob_mit != r.prob:
				ret += "\t\\begin{itemize}\n"
				if r.prob_mit != r.prob:
					ret += f"\t\t \item Probability mitigation (-{r.prob-r.prob_mit}): {r.prob_mit_exp.lower()}."
				if r.impact_mit != r.impact:
					ret += f"\t\t \item Impact mitigation (-{r.impact-r.impact_mit}): {r.impact_mit_exp.lower()}."
				ret += "\t\\end{itemize}\n"
		ret += "\end{itemize}"
		self.save_res(ret, "risk-def")
				
	def map(self):
		# 0=green, 1=yellow, 2=orange, 3=red
		color_map = [
			[3, 3, 3, 3, 3],
			[2, 2, 2, 3, 3],
			[0, 1, 1, 2, 3],
			[0, 0, 1, 2, 3],
			[0, 0, 0, 2, 3]
		]
		to_mitigate = 0
		tab_mitig = (", after mitigation", "-mitig") if self.use_mitig else ("", "")
		op = "\definecolor{rm-3}{HTML}{FE0000}\definecolor{rm-2}{HTML}{FF9900}\definecolor{rm-1}{HTML}{FCFF2F}\definecolor{rm-0}{HTML}{00FF00}"
		op += "\n\\begin{table}[H]\n"
		op += "\\centering\n"
		op += f"\\caption{{Risk map of the {self.ss_name} subsystem{tab_mitig[0]}}}\n"
		op += f"\\label{{tab:risk-map-{self.ss_codename.lower()}{tab_mitig[1]}}}"
		op += "\\begin{tabular}{l|c|c|c|c|c|}\n"
		op += "\\cline{2-6}\n"
		op += "& \multicolumn{1}{l|}{Very unlikely (1)} & \multicolumn{1}{l|}{Unlikely (2)} & \multicolumn{1}{l|}{Possible (3)} & \multicolumn{1}{l|}{Likely (4)} & \multicolumn{1}{l|}{Very likely (5)} \\\\ \hline" + "\n"
		l_impact = ["Very high", "High", "Medium", "Low", "Very Low"]
		for i in range(5):
			op += "\multicolumn{1}{|l|}{" + l_impact[i] + " impact (" + str(5-i) + ")}"
			t_r = self.get_by_impact(5-i)
			for j in range(5):
				in_cell = []
				for r in t_r:
					if self.use_mitig:
						if r.prob_mit == j+1:
							in_cell.append(str(r.id))
					else:
						if r.prob == j+1:
							in_cell.append(str(r.id))
				color = color_map[i][j]
				if color != 0 and len(in_cell) != 0:
					to_mitigate += 1
				op += f" & \cellcolor{{rm-{str(color)}}}" + ", ".join(in_cell)
			op += "\\\\ \\hline \n"
		op += "\\end{tabular} \n\\end{table}"
		self.save_res(op, f"risk-map{tab_mitig[1]}")
		if self.use_mitig and to_mitigate > 0:
			conjug = ("s are", "them") if to_mitigate > 1 else (" is", "it")
			print(f"/!\ {to_mitigate} risk{conjug[0]} still out of the green zone. Do your best to mitigate {conjug[1]} if possible /!\\")

	def save_res(self, res, fname, extension="txt"):
		f = open(self.path + f"\output\\{fname}.{extension}", "w")
		f.writelines(res)
		f.close()
		if not self.latex_file:
			print(f"Please copy the Latex risk map from output\\{fname}.txt to Overleaf.")
		else:
			self.latex.append(res)

	def mitig(self, status):
		self.use_mitig = status

	def combine_latex(self):
		out = f"\\noindent The following events have been assessed and mitigated as part of the {self.ss_name} subsystem:\n\n"
		out += self.latex[0]
		#out += f"\n\n\\noindent From this list, the risk map of \\autoref{{tab:risk-map-{self.ss_codename.lower()}}} has been created.\n\n"
		#out += self.latex[1]
		#out += f"\n\n\\noindent As seen in the risk map of \\autoref{{tab:risk-map-{self.ss_codename.lower()}}}, some risks have to be mitigated. \
		#An updated risk map, following mitigation, can be seen in \\autoref{{tab:risk-map-{self.ss_codename.lower()}-mitig}}.\n\n"
		out += f"\n\n\\noindent From this list, a mitigated risk map has been created, and can be seen in \\autoref{{tab:risk-map-{self.ss_codename.lower()}-mitig}}.\n\n"
		out += self.latex[2]
		out += "\n\n \\todo[inline]{Discuss the mitigated map, and alter the text if needed.}"
		self.save_res(out, f"risk-{self.ss_codename.lower()}", "tex")
		print(f"Please import the file output/risk-{self.ss_codename.lower()}.tex on Overleaf in the Chapters/Risks folder, and input it where needed.")