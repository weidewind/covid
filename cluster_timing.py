from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import datetime
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import numpy as np

parser = optparse.OptionParser()
parser.add_option('-a', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('-e', '--transmission_lineages_file', help='file with all strains assigned to lineage (ooutput from find_transmission_lineages.py)', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


def dateStats(cleaned_od):
	output_dates_formatted = np.array(pd.to_datetime(cleaned_od), dtype=np.datetime64)
	output_dates_num = pd.to_numeric(output_dates_formatted)

	median = pd.to_datetime(np.median(output_dates_num)).strftime("%Y-%m-%d")
	mean = pd.to_datetime(np.mean(output_dates_num)).strftime("%Y-%m-%d")
	min_date = pd.to_datetime(np.amin(output_dates_num)).strftime("%Y-%m-%d")
	max_date = pd.to_datetime(np.amax(output_dates_num)).strftime("%Y-%m-%d")

	return [mean, median, min_date, max_date]

print("Parsing dates")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))

with open(options.transmission_lineages_file, "r") as tlf:
	with open(options.output, "w") as out:
		out.write("\t".join(["entry", "mean_date", "median_date", "min_date", "min_date_strains","max_date", "all_dates"]) + "\n")
		for line in tlf:
			entry, strains = line.strip().split("\t")
			strain_list = strains.split(";")
			dates = [dates_dict[strain] for strain in strain_list if strain in dates_dict]
			cleaned_dates = [d for d in dates if d != "unknown" and len(d) >= 10]
			if cleaned_dates:
				fmean, fmedian, fmin, fmax = dateStats(cleaned_dates)
				earliest = [strain for strain in strain_list if dates_dict.get(strain,"unknown") == fmin]
				out.write("\t".join([entry, fmean, fmedian, fmin, ",".join(earliest), fmax, ",".join(cleaned_dates)]) + "\n")
			else:
				out.write("\t".join([entry, "", "", "", "", ""])+ "\n")