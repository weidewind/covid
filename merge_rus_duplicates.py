from meta import get_country_dict, get_gisaid_duplicates
import pandas as pd
import optparse

# add duplicates of russian strains from --duplicates to --rus_duplicates
# !!! unluckily, also add any nonrussian strains with russian duplicates 
# (=> the resulting file can contain nonrussian strains at the start of the line)

parser = optparse.OptionParser()
parser.add_option('-d', '--duplicates', help='file duplicates_21.03.txt, with true duplicates', type='str')
parser.add_option('-u', '--rus_duplicates', help='file all_duplicates_of_rus_strains, with kinda duplicates', type='str')
parser.add_option('-c', '--nonrus_clusters', help='file nonrus_clusters_far20.txt', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

meta_dict = get_country_dict()

dup_list = {}
centroid_for_strain = {}

print("Parsing --nonrus_clusters..")
with open(options.nonrus_clusters, "r") as clusts:
	for line in clusts:
		splitter = line.strip().split(';')
		# due to some bug there is one or several clusters containing russian strains
		if len(splitter) > 1 and "Russia" in set([meta_dict.get(s) for s in splitter[1:]]):
			rus_subset = [s for s in splitter[1:] if meta_dict.get(s) == "Russia"]
			dup_list[splitter[0]] = rus_subset
			for strain in rus_subset:
				centroid_for_strain[strain] = splitter[0]

print("Parsing --rus_duplicates..")
with open(options.rus_duplicates, "r") as rus_dupls:
	for line in rus_dupls:
		splitter = line.strip().split('\t')
		ds = splitter[1].split(';')
		if splitter[0] in dup_list:
			dup_list[splitter[0]].extend(ds)
		elif splitter[0] in centroid_for_strain:
			dup_list[centroid_for_strain[splitter[0]]].extend(ds)
			for strain in ds:
				centroid_for_strain[strain] = centroid_for_strain[splitter[0]]
		else:
			dup_list[splitter[0]] = ds


print("Parsing --duplicates..")
with open(options.duplicates, "r") as dupls:
	for line in dupls:
		splitter = line.strip().split(';')
		# there may be russian duplicates of nonrussian strains in this file. We also want them added into the tree
		if meta_dict.get(splitter[0]) == "Russia" or "Russia" in set([meta_dict.get(s) for s in splitter[1:]]):
			ds = splitter[1:]
			if splitter[0] in dup_list:
				dup_list[splitter[0]].extend(ds)
			elif splitter[0] in centroid_for_strain:
				dup_list[centroid_for_strain[splitter[0]]].extend(ds)
				for strain in ds:
					centroid_for_strain[strain] = centroid_for_strain[splitter[0]]
			else:
				dup_list[splitter[0]] = ds


gisaid_dupl = get_gisaid_duplicates()

for id in dup_list:
	if id in gisaid_dupl:
		print("Oh noes! Duplicate of russian strain, " + id + ", is one of the centroids!")
	elif id in meta['seq_id']:
		print("Oh noes! Russian strain with gisaid duplicate, " + id + ", is one of the centroids!")

print("Writing new duplicates list..")
with open(options.output, "w") as out:
	for k,v in dup_list.items():
		clean_v = [i for i in v if i not in gisaid_dupl]
		out.write(k + "\t" + ";".join(clean_v) + "\n")
