from meta import get_country_dict, get_gisaid_duplicates
import pandas as pd
import optparse

# mcollect duplicates of nonrus strains
# do not add nonrus duplicates that have rus duplicates

parser = optparse.OptionParser()
parser.add_option('-d', '--duplicates', help='file duplicates_21.03.txt, with true duplicates', type='str')
parser.add_option('-n', '--nonrus_duplicates', help='file nonrus.duplicates, with kinda duplicates', type='str')
parser.add_option('-c', '--nonrus_clusters', help='file nonrus_clusters_far20.txt', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

meta_dict = get_country_dict()

dup_list = {}
centroid_for_strain = {}

if (options.nonrus_clusters):
	print("Parsing --nonrus_clusters..")
	with open(options.nonrus_clusters, "r") as clusts:
		for line in clusts:
			splitter = line.strip().split(';')
			# due to some bug there is one or several clusters containing russian strains
			if len(splitter) > 1:
				nonrus_subset = [s for s in splitter[1:] if meta_dict.get(s) != "Russia"]
				dup_list[splitter[0]] = nonrus_subset
				for strain in nonrus_subset:
					centroid_for_strain[strain] = splitter[0]


print("Parsing --nonrus_duplicates..")

with open(options.nonrus_duplicates, "r") as nonrus_dupls:
	for line in nonrus_dupls:
		splitter = line.strip().split(';')
		if splitter[0] in dup_list:
			dup_list[splitter[0]].extend(splitter[1:])
			for strain in splitter[1:]:
				centroid_for_strain[strain] = splitter[0]
		elif splitter[0] in centroid_for_strain:
			dup_list[centroid_for_strain[splitter[0]]].extend(splitter[1:])
			for strain in splitter[1:]:
				centroid_for_strain[strain] = centroid_for_strain[splitter[0]]
		else:
			dup_list[splitter[0]] = splitter[1:]
			for strain in splitter[1:]:
				centroid_for_strain[strain] = splitter[0]


print("Parsing --duplicates..")
with open(options.duplicates, "r") as dupls:
	for line in dupls:
		splitter = line.strip().split(';')
		#if (splitter[0] == "EPI_ISL_428893"):
		#	print("I got him!\n")
		#	print("in meta_dict: " + meta_dict.get(splitter[0], "unknown"))
		# if there is a russian duplicate of non-russian strain, all these strains will be added while refluffing, and we don't want them here
		if not (meta_dict.get(splitter[0], "unknown") == "Russia" or "Russia" in set([meta_dict.get(s) for s in splitter[1:]])):
			if splitter[0] in dup_list:
				dup_list[splitter[0]].extend(splitter[1:])
				for strain in splitter[1:]:
					centroid_for_strain[strain] = splitter[0]
			elif splitter[0] in centroid_for_strain:
				dup_list[centroid_for_strain[splitter[0]]].extend(splitter[1:])
				for strain in splitter[1:]:
					centroid_for_strain[strain] = centroid_for_strain[splitter[0]]
			else:
				dup_list[splitter[0]] = splitter[1:]
				for strain in splitter[1:]:
					centroid_for_strain[strain] = splitter[0]


gisaid_dupl = get_gisaid_duplicates()

for id in dup_list:
	if id in gisaid_dupl:
		print("Duplicate of russian strain, " + id + ", is one of the centroids!")
	elif id in meta['seq_id']:
		print("Russian strain with gisaid duplicate, " + id + ", is one of the centroids!")

print("Writing new duplicates list..")
with open(options.output, "w") as out:
	for k,v in dup_list.items():
		clean_v = [i for i in v if i not in gisaid_dupl]
		out.write(k + "\t" + ";".join(clean_v) + "\n")
# gisaid_dupl = meta['gisaid_id'].tolist()
# with open (output, 'w') as out:
# 	for key, value in duplicates.items():
# 		clean_value = [i for i in value if i not in gisaid_dupl]
# 		out.write(key+"\t"+";".join(clean_value)+"\n")
