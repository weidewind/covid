from tree_utils import add_data_from_data_and_duplicates_files
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import re



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s', '--states', help='states file', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-d', '--duplicates', help='states file', type='str')
parser.add_option('-m', '--max_distance', help='maximum phylogenetic distance from entry node', type='float')
parser.add_option('-e', '--entry_nodes_file', help='file with entry nodes (one per line, output from find_transmission_lineages.py)', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()





def add_countries_from_countries_file(t, countries_dict):
	for node in t.traverse():
		if node.is_leaf() and re.match(r"^[a-zA-Z]", node.name):  # after pruning some internal nodes became leaves, hence the additional check
			if node.name in countries_dict:
				node.add_features(country=countries_dict[node.name])
			else:
				node.add_features(country="unknown")


def mark_and_collect_neighbours_countries(t, entry_node, node, dist, maxdist, countries_collection):
	if (not t.get_common_ancestor(node, entry_node) == node):
		dist += node.dist
	else:
		dist -= node.dist

	if (dist > maxdist):
		node.add_features(keep=False)
	else:
		node.add_features(keep=True)
		if not node == entry_node:
			if node.is_leaf():
				countries_collection.extend(node.country)
			else:
				for ch in node.children:
					mark_and_collect_neighbours_countries(t, entry_node, ch, dist, maxdist, countries_collection)


def mark_and_collect_neighbours(t, entry_node, node, dist, maxdist, neighbours_collection):
	if (not t.get_common_ancestor(node, entry_node) == node):
		dist += node.dist
	else:
		dist -= node.dist

	if (dist > maxdist):
		node.add_features(keep=False)
	else:
		node.add_features(keep=True)
		if not node == entry_node:
			if node.is_leaf():
				neighbours_collection.append(node.name)
			else:
				for ch in node.children:
					mark_and_collect_neighbours(t, entry_node, ch, dist, maxdist, neighbours_collection)


def mark_russian(t):
	for node in t.get_leaves():
		print("Marking " + node.name)
		if node.country == "Russia":
			node.add_features(keep=True)


def multiple_mark_russian(t):
	for node in t.get_leaves():
		print("Marking " + node.name)
		if "Russia" in node.country:
			node.add_features(keep=True)


def strip_tree(t, outfile):
	# left, then right, then root
	for node in t.traverse("postorder"):
		if hasattr(node, 'keep') and node.keep == True:
			if not node.is_root():
				node.up.add_features(keep=True)
				for ch in node.up.get_children():
					ch.add_features(keep=True)
		else:
			node.add_features(keep=False)

	t.write(is_leaf_fn=last_kept, format=1, outfile=outfile, format_root_node=True)
	t2 = Tree(outfile, format=1)
	return (t2)


# should this node be considered a leaf?
def last_kept(node):
	if node.is_leaf():
		return(True)
	elif all(not ch.keep for ch in node.get_children()):
		return(True)
	else:
		return(False)


def color_tree(t):
	for node in t.traverse():
		print(node.name)
		nstyle = NodeStyle()
		nstyle["fgcolor"] = matplotlib.colors.to_hex([states[node.name], 0.0, 0.0])
		nstyle["size"] = 1
		if states[node.name] > 0:
			nstyle["size"] = 2
		if node.name in entries:
			nstyle["size"] = 3
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		node.set_style(nstyle)


def countries_collection_to_str(c_collection):
	stats = dict(Counter(c_collection))
	sorted_stats = sorted(stats.items(), key=lambda x: x[1], reverse=True)
	string = ";".join([c[0] + " " + str(c[1]) for c in sorted_stats])
	return (string)


def country_layout(node):
	if node.is_leaf() and re.match(r"^[a-zA-Z]", node.name) and not node.country == "Russia":
		faces.add_face_to_node(AttrFace("country"), node, column=0)
	if node.name in entries:
		faces.add_face_to_node(AttrFace("name"), node, column=0)


def multiple_country_layout(node):
	if node.is_leaf() and re.match(r"^[a-zA-Z]", node.name):
		faces.add_face_to_node(faces.TextFace(countries_collection_to_str(node.country)), node, column=0)
	if node.name in entries:
		faces.add_face_to_node(AttrFace("name"), node, column=0)


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

print("Assigning locations..")
add_data_from_data_and_duplicates_files(tree, countries_dict, options.duplicates, "country")
print("Marking russian strains so they won't be pruned..")
multiple_mark_russian(tree)

print("Collecting neighbours.. ")
entries = []
with open(options.entry_nodes_file, "r") as enf:
	with open(options.output + ".country_stats", "w") as out:
		out.write("\t".join(["entry", "closest_countries_stats", "closest_countries", "closest_nodes"]) + "\n")
		for line in enf:
			print(">entry_node " + line.strip())
			out.write(line.strip() + "\t")
			entry_node = tree.search_nodes(name=line.strip())[0]
			entries.append(entry_node.name)
			node = entry_node
			dist = node.dist
			while dist < options.max_distance:
				node = node.up
				dist += node.dist #its ok to add the last distance too - we will subtract it later, in mark_and_collect_neighbours
			neighbours_collection = []
			mark_and_collect_neighbours(t=tree, entry_node=entry_node, node=node, dist=dist, maxdist=options.max_distance, neighbours_collection=neighbours_collection)
			strains_collection = []
			for n in neighbours_collection:
				strains_collection.append(n)
				if n in duplicates:
					for d in duplicates[n]:
						strains_collection.append(d)
			countries_collection = [countries_dict[n] for n in neighbours_collection]
			stats = dict(Counter(countries_collection))
			stats_str = ",".join([k + ":" + str(v) for k,v in stats.items()])
			out.write("\t".join([stats_str, ",".join(countries_collection), ",".join([n for n in neighbours_collection])]) + "\n")

print("Stripping tree..")
t2 = strip_tree(tree, options.output)
t2.write(format=1, outfile=options.output + "_test2")

print("Parsing states..")
states = {}
with open(options.states) as st:
	for line in st:
		[nname, _, state2_prob] = line.split()
		states[nname] = float(state2_prob)  # probability of being russian
		print("nname " + nname + " state 2 probability " + str(state2_prob))

print("Styling tree..")
add_data_from_data_and_duplicates_files(t2, countries_dict, options.duplicates, "country")
color_tree(t2)

ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = multiple_country_layout
ts.scale = 1
t2.ladderize()
#t2.render(options.output + ".png", tree_style=ts, dpi=300, w=870, units="mm")
t2.render(options.output + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
t2.render(options.output + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")


