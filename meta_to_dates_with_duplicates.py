from meta import get_date_dict
from ete3 import Tree
import pandas as pd
import numpy as np
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--duplicates', help='file with duplicates (and cluster memebers)', type='str')
parser.add_option('-o', '--output', help='output', type='str')
parser.add_option('-y', '--type', help='min, max, mean or median of all dates?', type='str')

options, args = parser.parse_args()


def dateStats(cleaned_od, quanttype):
	if len(cleaned_od) == 0:
		return ""
	
	output_dates_formatted = np.array(pd.to_datetime(cleaned_od, errors='coerce'), dtype=np.datetime64)
	output_dates_formatted = output_dates_formatted[~ np.isnat(output_dates_formatted)]
	output_dates_num = pd.to_numeric(output_dates_formatted)

	if quanttype == "median":
		return pd.to_datetime(np.median(output_dates_num)).strftime("%Y-%m-%d")
	if quanttype == "mean":
		return pd.to_datetime(np.mean(output_dates_num)).strftime("%Y-%m-%d")
	if quanttype == "min":
		try:
				print(pd.to_datetime(np.amin(output_dates_num)))
		except:
				print("Cleaned_od: ")
				print(",".join(cleaned_od))
				print("output_dates_formatted:")
				print(",".join(output_dates_formatted))
				print("output_dates_num:")
				print(",".join(output_dates_num))
				raise ValueError('oops')
		return pd.to_datetime(np.amin(output_dates_num)).strftime("%Y-%m-%d")
	if quanttype == "max":
		return pd.to_datetime(np.amax(output_dates_num)).strftime("%Y-%m-%d")
	else:
		raise ValueError('Undefined type: ' + quanttype)


print("Parsing tree..")
tree = Tree(options.tree, format=1)

meta_dict = get_date_dict()


print("Parsing duplicates..")
duplicates = {}
with open(options.duplicates, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if len(splitter) > 1:
			duplicates[splitter[0]] = splitter[1].split(';')


print("Writing output..")
with open(options.output + "." + options.type, "w") as out:
	out.write("name,date" + "\n")
	for n in tree.iter_leaves():
		strain = n.name
		date = meta_dict.get(strain, "")
		if options.type == "all":
			out.write(strain + "," + date + "\n")
			if strain in duplicates:
				for dup in duplicates[strain]:
					out.write(dup + "," + meta_dict.get(dup, "") + "\n")
		else:
			s_dates = []
			s_dates.append(date)
			if strain in duplicates:
				for dup in duplicates[strain]:
					s_dates.append(meta_dict.get(dup, "unknown"))
			cleaned_s_dates = [d for d in s_dates if d != "unknown" and len(d) == 10 ] # exclude dates with unknown exact day of month
			out.write(strain + "," + dateStats(cleaned_s_dates, options.type) + "\n")
