import qualitative as ql

"""
How to use:
	1: Change the name of your susbsystem (ss_name), as well as its ~4 letters codename (ss_codename)
	2: Duplicate the Example-COPY-ONLY csv file in the input folder, change its name, and update it to match your own risk analysis
	3: Run this file, import your output "tex" file on Overleaf, and \include{Chapters/risk/your-risk.tex}
	4: Adapt the Overleaf "tex" file to match the maps, and discuss the results
"""

files = {
	"prop": ["Propulsion", "propulsion-risk"],
	"thermal": ["Thermal", "thermal_risk"],
	"gnc": ["GNC", "GNC_risks"],
	"struc": ["Structures", "Structures-risk"],
	"power": ["Power", "Power-risk"]
	}

for key, value in files.items():
	risks = ql.risks_list(ss_codename=key, ss_name=value[0], risk_list=[])	# Input subsystem codename
	risks.get_risks(fname=value[1])	# Change input file (copy example and change the values)