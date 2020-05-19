import qualitative as ql

"""
How to use:
	1: Change the name of your susbsystem (ss_name), as well as its ~4 letters codename (ss_codename)
	2: Duplicate the Example-COPY-ONLY csv file in the input folder, change its name, and update it to match your own risk analysis
	3: Run this file, import your output "tex" file on Overleaf, and \include{Chapters/risk/your-risk.tex}
	4: Adapt the Overleaf "tex" file to match the maps, and discuss the results
"""
risks = ql.risks_list(ss_codename="safe", ss_name="Life Support")		# Input subsystem codename
risks.get_risks(fname="Example-COPY-ONLY")	# Change input file (copy example and change the values)