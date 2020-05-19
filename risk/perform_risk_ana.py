import qualitative as ql

risks = ql.risks_list(ss_codename="safe", ss_name="Life Support")		# Input subsystem codename
risks.get_risks(fname="Example-COPY-ONLY")	# Change input file (copy example and change the values)