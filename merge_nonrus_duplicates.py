from meta import get_country_dict, get_gisaid_duplicates
import pandas as pd
import optparse

parser = optparse.OptionParser()
parser.add_option('-n', '--nonrus_duplicates', help='file nonrus.duplicates, with kinda duplicates', type='str')
parser.add_option('-c', '--nonrus_clusters', help='file nonrus_clusters_far20.txt', type='str')
parser.add_option('-w', '--without_clusters', action="store_true", default=False)
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

meta_dict = get_country_dict()

dup_list = {}
centroid_for_strain = {}

if not options.without_clusters:
	print("Parsing --nonrus_clusters..")
	with open(options.nonrus_clusters, "r") as clusts:
		for line in clusts:
			splitter = line.strip().split('\t')
			if len(splitter) > 1:
				strains = splitter[1].split(',')
				dup_list[splitter[0]] = strains
				for strain in strains:
					centroid_for_strain[strain] = splitter[0]


print("Parsing --nonrus_duplicates..")

with open(options.nonrus_duplicates, "r") as nonrus_dupls:
	for line in nonrus_dupls:
		splitter = line.strip().split('\t')
		lead = splitter[0]
		dups = splitter[1].split(";")
		if lead in dup_list:
			dup_list[lead].extend(dups)
			for strain in dups:
				centroid_for_strain[strain] = lead
		elif lead in centroid_for_strain:
			dup_list[centroid_for_strain[lead]].extend(dups)
			for strain in dups:
				centroid_for_strain[strain] = centroid_for_strain[lead]
		else:
			dup_list[lead] = dups
			for strain in dups:
				centroid_for_strain[strain] = lead


print("Writing new duplicates list..")
with open(options.output, "w") as out:
	for k,v in dup_list.items():
		out.write(k + "\t" + ";".join(v) + "\n")
# gisaid_dupl = meta['gisaid_id'].tolist()
# with open (output, 'w') as out:
# 	for key, value in duplicates.items():
# 		clean_value = [i for i in value if i not in gisaid_dupl]
# 		out.write(key+"\t"+";".join(clean_value)+"\n")
