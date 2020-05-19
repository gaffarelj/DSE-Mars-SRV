import qualitative as ql

risks = ql.risks_list(ss_name="safe")		# Input subsystem codename
risks.get_risks(fname="Example-COPY-ONLY")	# Change input file (copy example and change the values)
risks.risks_def()							# Define the risks before the maps
risks.map()									# Generate a map before mitigation
risks.mitig(True)
risks.map()									# Generate a map after mitigation