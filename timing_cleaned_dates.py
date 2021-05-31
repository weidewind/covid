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
parser.add_option('-e', '--entry_nodes_file', help='file with entry nodes (one per line, output from find_transmission_lineages.py)', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


# def leaf_not_fully_russian(leaf):
# 	country_set = list(set(leaf.country))
# 	if len(country_set) > 1 or not country_set[0] == "Russia":
# 		return True
# 	else:
# 		return False

def leaf_not_fully_russian(leaf):
	if countries_dict[leaf.name] not in ["Russia", "unknown"]:
		return True
	if leaf.name in duplicates:
		for s in duplicates[leaf.name]:
			if s in countries_dict and countries_dict[s] not in ["Russia", "unknown"]:
				return True
		return False


def leaf_has_russians(leaf):
	if countries_dict[leaf.name] == "Russia":
		return True
	if leaf.name in duplicates:
		for s in duplicates[leaf.name]:
			if s in countries_dict and countries_dict[s] == "Russia":
				return True
		return False


def collectClosest(node, curdist, mindist, output_nodes, leaf_is_ok):  # curdepth - current distance to the leaf in question; mindist - best distance so far
	print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_foreign_nodes", ",".join([n.name for n in output_nodes])]))
	if node.is_leaf() and curdist <= mindist and leaf_is_ok(node):
		if curdist < mindist:
			output_nodes.clear()
		output_nodes.append(node)
		mindist = curdist
	elif curdist <= mindist:
		for ch in node.children:
			chdist = curdist + ch.dist
			if chdist <= mindist:
				(mindist, output_nodes) = collectClosest(ch, chdist, mindist, output_nodes, leaf_is_ok)
		print("After traversing all children of node " + node.name + ":")
		print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_foreign_nodes", ",".join([n.name for n in output_nodes])]))
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


def collectAndCleanDates(nodes, country_is_ok):
	output_dates = []
	for n in nodes:
		if country_is_ok(countries_dict[n.name]):
			output_dates.append(dates_dict[n.name])
		if n.name in duplicates:
			for d in duplicates[n.name]:
				if d in countries_dict and country_is_ok(countries_dict[d]):
					output_dates.append(dates_dict[d])
	cleaned_od = [d for d in output_dates if d != "unknown"]
	return cleaned_od, len(output_dates)


def dateStats(cleaned_od):
	output_dates_formatted = np.array(pd.to_datetime(cleaned_od), dtype=np.datetime64)
	output_dates_num = pd.to_numeric(output_dates_formatted)

	median = pd.to_datetime(np.median(output_dates_num)).strftime("%Y-%m-%d")
	mean = pd.to_datetime(np.mean(output_dates_num)).strftime("%Y-%m-%d")
	min_date = pd.to_datetime(np.amin(output_dates_num)).strftime("%Y-%m-%d")
	max_date = pd.to_datetime(np.amax(output_dates_num)).strftime("%Y-%m-%d")

	return [mean, median, min_date, max_date]


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

# print("Assigning dates..")
# add_data_from_data_and_duplicates_files(tree, dates_dict, options.duplicates, "date")
# print("Assigning countries..")
# add_data_from_data_and_duplicates_files(tree, countries_dict, options.duplicates, "country")

with open(options.entry_nodes_file, "r") as enf:
	with open(options.output + ".dates_stats", "w") as out:
		out.write("\t".join(["entry", "mean_foreign_date", "median_foreign_date", "min_foreign_date", "max_foreign_date", "number_of_close_foreign_strains", "close_foreign_dates", "close_foreign_strains", "tree_distance_to_closest_foreign_strain"]) + "\t")
		out.write("\t".join(["mean_rus_date", "median_rus_date", "min_rus_date", "max_rus_date", "number_of_close_rus_strains", "close_rus_dates", "close_rus_strains", "tree_distance_to_closest_rus_strain"]) + "\n")

		for line in enf:
			print(">entry_node " + line.strip())
			out.write(line.strip() + "\t")
			entry_node = tree.search_nodes(name=line.strip())[0]
			farthest, maxdist = tree.get_farthest_node()

			#  process closest foreign nodes
			output_foreign_nodes = []
			tnode = entry_node
			fmindist = maxdist

			# add the entry node itself, since it could also have some foreign duplicates
			if entry_node.is_leaf() and leaf_not_fully_russian(entry_node):
				output_foreign_nodes.append(entry_node)
				fmindist = 0
			else:
				for ch in entry_node.children:
					if ch.is_leaf() and ch.dist == 0 and leaf_not_fully_russian(ch):
						output_foreign_nodes.append(ch)	
						fmindist = 0

			while(not tnode.is_root() and entry_node.get_distance(tnode) <= fmindist):
				siss = tnode.get_sisters()
				for s in siss:
					print("sis " + s.name)
					curdist = entry_node.get_distance(s)
					fmindist, output_foreign_nodes = collectClosest(s, curdist, fmindist, output_foreign_nodes, leaf_not_fully_russian)
				tnode = tnode.up

			print("Collected foreign nodes:")
			print(",".join([n.name for n in output_foreign_nodes]))
			

			cleaned_foreign_dates, dirty_foreign_num = collectAndCleanDates(output_foreign_nodes, country_is_foreign)
			if cleaned_foreign_dates:
				fmean, fmedian, fmin, fmax = dateStats(cleaned_foreign_dates)

			# process closest russian nodes
			output_rus_nodes = []
			rmindist = maxdist
			rmindist, output_rus_nodes = collectClosest(entry_node, 0, rmindist, output_rus_nodes, leaf_has_russians)

			print("Collected rus nodes:")
			print(",".join([n.name for n in output_rus_nodes]))

			cleaned_rus_dates, dirty_rus_num = collectAndCleanDates(output_rus_nodes, country_is_russia)
			if cleaned_rus_dates:
				rmean, rmedian, rmin, rmax = dateStats(cleaned_rus_dates)

			out.write("\t".join([fmean, fmedian, fmin, fmax, str(dirty_foreign_num), ",".join(cleaned_foreign_dates), ",".join([n.name for n in output_foreign_nodes]), str(fmindist)])+ "\t")
			out.write("\t".join([rmean, rmedian, rmin, rmax, str(dirty_rus_num), ",".join(cleaned_rus_dates), ",".join([n.name for n in output_rus_nodes]), str(rmindist)]) + "\n")
