from tree_utils import add_data_from_data_and_duplicates_files
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import datetime
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import numpy as np
import re

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-a', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-d', '--duplicates', help='duplicates file', type='str')
parser.add_option('-m', '--maxdist', help='maximum distance for closest nodes', type='float', default = None)
parser.add_option('-e', '--entry_nodes_file', help='file with entry nodes (one per line, output from find_transmission_lineages.py)', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

def get_distance_to_parent(child, parent):
	temp = child
	dist = 0
	while not temp.name == parent.name:
		dist += temp.dist
		temp = temp.up
	return(dist)
# def leaf_not_fully_russian(leaf):
# 	country_set = list(set(leaf.country))
# 	if len(country_set) > 1 or not country_set[0] == "Russia":
# 		return True
# 	else:
# 		return False

# def leaf_not_fully_russian(leaf):
# 	if countries_dict[leaf.name] not in ["Russia", "unknown"]:
# 		return True
# 	if leaf.name in duplicates:
# 		for s in duplicates[leaf.name]:
# 			if s in countries_dict and countries_dict[s] not in ["Russia", "unknown"]:
# 				return True
# 		return False


# def leaf_has_russians(leaf):
# 	if countries_dict[leaf.name] == "Russia":
# 		return True
# 	if leaf.name in duplicates:
# 		for s in duplicates[leaf.name]:
# 			if s in countries_dict and countries_dict[s] == "Russia":
# 				return True
# 		return False


def leaf_has_dated_russian(leaf):
	if len(node_to_dates[leaf.name]["r"]) > 0:
		return True
	else:
		return False


def leaf_has_dated_foreign(leaf):
	if len(node_to_dates[leaf.name]["f"]) > 0:
		return True
	else:
		return False


def any_leaf(leaf):
	return True


def collectClosest(node, curdist, mindist, output_nodes, leaf_is_ok):  # curdepth - current distance to the leaf in question; mindist - best distance so far
	#print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_nodes", ",".join([n.name for n in output_nodes])]))
	if node.is_leaf() and curdist <= mindist and leaf_is_ok(node):
		if curdist < mindist:
			output_nodes.clear()
		output_nodes.append(node)
		mindist = curdist
	elif curdist <= mindist:
		for ch in node.children:
			chdist = curdist + ch.dist
			if chdist <= mindist: ## add or output_nodes does not contain dates
				(mindist, output_nodes) = collectClosest(ch, chdist, mindist, output_nodes, leaf_is_ok)
		#print("After traversing all children of node " + node.name + ":")
		#print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_nodes", ",".join([n.name for n in output_nodes])]))
	return(mindist, output_nodes)


def country_is_russia(country):
	if country == "Russia":
		return True
	else:
		return False


def country_is_foreign(country):
	if country not in ["Russia", "unknown"]:
		return True
	else:
		return False

def any_country(country):
	return True


# def collectDates(nodes, country_is_ok):
# 	output_dates = []
# 	corresponding_nodes = []
# 	add_nodes_date = [dates_dict[n.name] for n in nodes if country_is_ok(countries_dict[n.name])]
# 	output_dates.extend(add_nodes_date)
# 	corresponding_nodes.extend([n for n in nodes if country_is_ok(countries_dict[n.name])])
# 	for n in nodes:
# 		if n.name in duplicates:
# 			add_dates = [dates_dict.get(d,"unknown") for d in duplicates[n.name] if country_is_ok(countries_dict.get(d, "unknown"))]
# 			output_dates.extend(add_dates)
# 			corresponding_nodes.extend([n] * len(add_dates))
# 	return output_dates, corresponding_nodes


def collectDates(nodes, country_pointer):
	output_dates = []
	corresponding_nodes = []
	for n in nodes:
		dates = node_to_dates[n.name][country_pointer]
		output_dates.extend(dates)
		corresponding_nodes.extend([n] * len(dates))
	return output_dates, corresponding_nodes


def dateStats(cleaned_od):
	output_dates_formatted = np.array(pd.to_datetime(cleaned_od), dtype=np.datetime64)
	output_dates_num = pd.to_numeric(output_dates_formatted)

	median = pd.to_datetime(np.median(output_dates_num)).strftime("%Y-%m-%d")
	mean = pd.to_datetime(np.mean(output_dates_num)).strftime("%Y-%m-%d")
	min_date = pd.to_datetime(np.amin(output_dates_num)).strftime("%Y-%m-%d")
	max_date = pd.to_datetime(np.amax(output_dates_num)).strftime("%Y-%m-%d")

	return [mean, median, min_date, max_date]


def collectClosestSearch(entry_node, leaf_is_ok, country_pointer, maxdist=None):
	output_nodes = []
	if maxdist is None:
		farthest, maxdist = tree.get_farthest_node()
	mindist = maxdist
	tnode = entry_node

	mindist, output_nodes = collectClosest(entry_node, 0, mindist, output_nodes, leaf_is_ok)

	while(not tnode.is_root() and get_distance_to_parent(child=entry_node, parent=tnode) <= mindist):
		siss = tnode.get_sisters()
		for s in siss:
			#print("sis " + s.name)
			curdist = get_distance_to_parent(child=entry_node, parent=tnode) + tnode.dist + s.dist
			mindist, output_nodes = collectClosest(s, curdist, mindist, output_nodes, leaf_is_ok)
		tnode = tnode.up

	print("Finished.")
	print("Collecting dates and distances..")
	#print(",".join([n.name for n in output_nodes]))
	clean_dates, corresponding_clean_nodes = collectDates(output_nodes, country_pointer)
	return clean_dates, corresponding_clean_nodes, mindist


def get_mrca_distances(entry_node, nodes, entry_dist_dict):
	dist_mrca = []
	dist_mrca_entry = []
	entry_dist_dict = {}
	#print("Dirty outgroup dates " + ",".join(dirty_foreign_dates))
	#print("Corr dirty outgroup nodes " + ",".join([fn.name for fn in corresponding_dirty_foreign_nodes]))
	for n in nodes:
		if n.name not in entry_dist_dict:
			mrca = tree.get_common_ancestor(n, entry_node)
			dm = get_distance_to_parent(child = n, parent = mrca)
			dme = get_distance_to_parent(child = entry_node, parent = mrca)
			dist_mrca.append(dm)
			dist_mrca_entry.append(dme)
			entry_dist_dict[n.name] = [dm, dme]
		else:
			dm, dme = entry_dist_dict[n.name]
			dist_mrca.append(dm)
			dist_mrca_entry.append(dme)
	return dist_mrca_entry, dist_mrca


print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing dates")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))

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

print("Connecting nodes, countries and dates..")
node_to_dates = {} # node_to_dates{leaf_name}{"r"}["2020-03-10", "2020-04-11"] - all correct dates for russia(r) and other countries(f)
for leaf in tree.iter_leaves():
	fstrains = []
	rstrains = []
	node_to_dates[leaf.name] = {"f":[], "r":[]}
	strains = [leaf.name]
	if leaf.name in duplicates:
		strains.extend(duplicates[leaf.name])
	for s in strains:
		c = countries_dict.get(s, "unknown")
		d = dates_dict.get(s, "unknown")
		if d != "unknown" :
			if c == "Russia": 
				node_to_dates[leaf.name]["r"].append(d)
			else:
				node_to_dates[leaf.name]["f"].append(d)


print("Processing entries..")
with open(options.entry_nodes_file, "r") as enf:
	with open(options.output + ".dates_stats", "w") as out:
		out.write("\t".join(["entry", "mean_foreign_date", "median_foreign_date", "min_foreign_date", "max_foreign_date", "number_of_close_foreign_strains", "close_foreign_dates", "close_foreign_strains", "tree_distance_foreign_to_mrca", "tree_distance_foreign_mrca_to_entry", "tree_distance_to_closest_foreign_strain"]) + "\t")
		out.write("\t".join(["mean_rus_date", "median_rus_date", "min_rus_date", "max_rus_date", "number_of_close_rus_strains", "close_rus_dates", "close_rus_strains", "tree_distance_rus_to_mrca", "tree_distance_rus_mrca_to_entry", "tree_distance_to_closest_rus_strain"]) + "\n")

		for line in enf:
			print("Collecting closest nodes with dated foreign strains..")
			print(">entry_node " + line.strip())
			out.write(line.strip() + "\t")
			entry_node = tree.search_nodes(name=line.strip())[0]
			entry_dist_dict = {}
			
			foreign_clean_dates, foreign_clean_nodes, fmindist = collectClosestSearch(entry_node, leaf_has_dated_foreign, country_pointer="f", maxdist=options.maxdist)
			foreign_dist_mrca_entry, foreign_dist_mrca = get_mrca_distances(entry_node, foreign_clean_nodes, entry_dist_dict)
			print("entry_dist_dict:")
			print(",".join([k + str(v) for k,v in entry_dist_dict.items()]))
			foreign_clean_num = len(foreign_clean_dates)
			if foreign_clean_dates:
				fmean, fmedian, fmin, fmax = dateStats(foreign_clean_dates)
			else:
				print("No foreign_clean_dates for " + entry_node.name + "! Will stop here.")
				raise

			# process closest russian nodes
			print("Collecting closest nodes with dated rus strains..")

			rus_clean_dates, rus_clean_nodes, rmindist = collectClosestSearch(entry_node, leaf_has_dated_russian, country_pointer="r", maxdist=options.maxdist)
			rus_dist_mrca_entry, rus_dist_mrca = get_mrca_distances(entry_node, rus_clean_nodes, entry_dist_dict)

			rus_clean_num = len(rus_clean_dates)
			if rus_clean_dates:
				rmean, rmedian, rmin, rmax = dateStats(rus_clean_dates)
			else:
				print("No clean_rus_dates for " + entry_node.name + "! Will stop here.")
				raise

			out.write("\t".join([fmean, fmedian, fmin, fmax, str(foreign_clean_num), ",".join(foreign_clean_dates), ",".join([n.name for n in foreign_clean_nodes]), ",".join([str(d) for d in foreign_dist_mrca]), ",".join([str(d) for d in foreign_dist_mrca_entry]), str(fmindist)])+ "\t")
			out.write("\t".join([rmean, rmedian, rmin, rmax, str(rus_clean_num), ",".join(rus_clean_dates), ",".join([n.name for n in rus_clean_nodes]), ",".join([str(d) for d in rus_dist_mrca]), ",".join([str(d) for d in rus_dist_mrca_entry]), str(rmindist)]) + "\n")
