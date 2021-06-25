from tree_utils import add_data_from_data_and_duplicates_files
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import re

# python country_stat.py --duplicates ../data/munched/gennady/may31/21_05_2021_all_nonrus_duplicates 
#--tree ../output/transmission_lineages/gennady/may31/Sankoff_states.csv.newick 
#--states ../output/transmission_lineages/gennady/may31/Sankoff_states.csv.probs 
#--countries ../output/transmission_lineages/gennady/may31/leaf_states.csv
# --entry_nodes_file ../output/transmission_lineages/gennady/may31/transmission_lineages.withduplicates.out.entries 
#--countries ../output/transmission_lineages/gennady/may31/leaf_states.csv
# --max_distance 0.00012
# --output ../output/transmission_lineages/gennady/may31/stripped_tree_with_labels_with_duplicates.00012
# >../output/transmission_lineages/gennady/may31/stripped_tree_with_labels_with_duplicates.logger.00012

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s', '--states', help='states file', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-d', '--duplicates', help='states file', type='str')
parser.add_option('-m', '--max_distance', help='maximum phylogenetic distance from entry node', type='float')
parser.add_option('-e', '--entry_nodes_file', help='file with entry nodes (one per line, output from find_transmission_lineages.py)', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


def collect_neighbours(t, entry_node, node, dist, maxdist, neighbours_collection, distances_collection, mrca_collection, entry_parents, mrca):
	if (node.name not in entry_parents):
		dist += node.dist
	else:
		dist -= node.dist
		mrca = node

	if (dist <= maxdist):
		if not node == entry_node:
			if node.is_leaf():
				neighbours_collection.append(node.name)
				distances_collection.append(dist)
				mrca_collection.append(mrca.name)
				#print("Appended " + node.name + " at " + str(dist) + " with mrca " + mrca.name)

			else:
				for ch in node.children:
					collect_neighbours(t, entry_node, ch, dist, maxdist, neighbours_collection, distances_collection, mrca_collection, entry_parents, mrca)



print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing countries..")
countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
countries_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(countries_dict.items())[0:5]))

print("Parsing duplicates..")
duplicates = {}
with open(options.duplicates, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if len(splitter) > 1:
			duplicates[splitter[0]] = splitter[1].split(';')


print("Collecting neighbours.. ")
entries = []
with open(options.entry_nodes_file, "r") as enf:
	with open(options.output + ".country_stats", "w") as out:
		out.write("\t".join(["entry", "number_of_closest_strains", "closest_countries_stats", "closest_countries", "closest_nodes", "distances", "etm_distances"]) + "\n")
		for line in enf:
			print(">entry_node " + line.strip())
			out.write(line.strip() + "\t")
			entry_node = tree.search_nodes(name=line.strip())[0]
			entries.append(entry_node.name)
			node = entry_node
			dist = node.dist
			mrca_dist_dict = {}
			while dist < options.max_distance and not node.is_root():
				node = node.up
				mrca_dist_dict[node.name] = dist
				dist += node.dist #its ok to add the last distance too - we will subtract it later, in mark_and_collect_neighbours
				

			entry_parents = [entry_node.name]
			pnode = entry_node
			while not pnode.is_root():
				pnode = pnode.up
				entry_parents.append(pnode.name)


			neighbours_collection = []
			distances_collection = []
			mrca_collection = []
			collect_neighbours(t=tree, entry_node=entry_node, node=node, dist=dist, maxdist=options.max_distance, distances_collection=distances_collection, neighbours_collection=neighbours_collection, mrca_collection=mrca_collection, mrca=node, entry_parents=set(entry_parents))
			strains_collection = []
			strains_dist_collection = []
			strains_etmdist_collection = []
			strains_collection.extend(neighbours_collection)
			strains_dist_collection.extend(distances_collection)
			strains_etmdist_collection.extend([mrca_dist_dict[mname] for mname in mrca_collection])
			
			for ind, n in enumerate(neighbours_collection):
				if n in duplicates:
					dups = duplicates[n]
					strains_collection.extend(dups)
					print("dups " + ",".join(dups))
					print("ind " + str(ind))
					print("dist " + str(distances_collection[ind]))
					strains_dist_collection.extend([distances_collection[ind]]*len(dups))
					strains_etmdist_collection.extend([mrca_dist_dict[mrca_collection[ind]]]*len(dups))

			countries_collection = [countries_dict.get(n, "unknown") for n in strains_collection]
			stats = dict(Counter(countries_collection))
			stats_str = ",".join([k + ":" + str(v) for k,v in stats.items()])
			out.write("\t".join([str(len(countries_collection)), stats_str, ",".join(countries_collection), ",".join([n for n in strains_collection]), ",".join([str(round(n,11)) for n in strains_dist_collection]), ",".join([str(round(n,11)) for n in strains_etmdist_collection])]) + "\n")

