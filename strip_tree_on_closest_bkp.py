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
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-a', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('-d', '--duplicates', help='states file', type='str')
parser.add_option('-e', '--entry_nodes_file', help='file with entry nodes (one per line, output from find_transmission_lineages.py)', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


def leaf_not_fully_russian(leaf):
	if countries_dict[leaf.name] not in ["Russia", "unknown"]:
		return True
	if leaf.name in duplicates:
		for s in duplicates[leaf.name]:
			if s in countries_dict and countries_dict[s] not in ["Russia", "unknown"]:
				return True
		return False


def get_distance_to_parent(child, parent):
	temp = child
	dist = 0
	while not temp.name == parent.name:
		dist += temp.dist
		temp = temp.up
	return(dist)


def collectClosest(node, curdist, mindist, output_nodes, leaf_is_ok):  # curdepth - current distance to the leaf in question; mindist - best distance so far
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
	return(mindist, output_nodes)


def strip_tree(t, outfile):
	# left, then right, then root
	for node in t.traverse("postorder"):
		if hasattr(node, 'keep') and node.keep:
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


def countries_collection_to_str(c_collection):
	stats = dict(Counter(c_collection))
	sorted_stats = sorted(stats.items(), key=lambda x: x[1], reverse=True)
	string = ";".join([c[0] + " " + str(c[1]) for c in sorted_stats])
	return (string)


def multiple_country_layout(node):
	if node.is_leaf() and re.match(r"^[a-zA-Z]", node.name) and node.name not in entries:
		faces.add_face_to_node(faces.TextFace(countries_collection_to_str(node.country)), node, column=0)
	if node.name in entries:
		faces.add_face_to_node(faces.TextFace(node.name + " " + str(node.countries_str), fgcolor="#ee0000"), node,  column=0)
		if node.count > 1:
			faces.add_face_to_node(faces.TextFace("earliest " + node.min_date[:11], fgcolor="#ee0000"), node,  column=0)
			faces.add_face_to_node(faces.TextFace("latest " + node.max_date[:11], fgcolor="#ee0000"), node,  column=0)
		else:
			faces.add_face_to_node(faces.TextFace(node.min_date[:11], fgcolor="#ee0000"), node,  column=0)
		#node.add_face(faces.TextFace(node.name + " " + str(node.countries_str)),column=0)
		#faces.add_face_to_node(AttrFace("name"), node, column=0)
		#faces.add_face_to_node(AttrFace("strains_str"), node, column=1)
		#faces.add_face_to_node(AttrFace("countries_str"), node, column=2)
		#print(node.strains_str)
		#node.add_face(faces.TextFace(node.name + " " + str(node.strains_str)),column=0) # ok but doubled
		#node.add_face(faces.TextFace(node.strains_str),column=1)
		#faces.add_face_to_node(faces.TextFace(node.name), node, column=0)
		#faces.add_face_to_node(faces.TextFace(node.strains_str), node, column=0)
		#faces.add_face_to_node(faces.TextFace(node.countries_str), node, column=2)


def color_tree(t):
	for node in t.traverse():
		print(node.name)
		nstyle = NodeStyle()
		nstyle["size"] = 1
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		if node.name in entries:
			nstyle["fgcolor"] = "#ee0000"
			nstyle["size"] = 4
		node.set_style(nstyle)

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing countries..")
countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
countries_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(countries_dict.items())[0:5]))

print("Parsing dates..")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))

print("Parsing duplicates..")
duplicates = {}
with open(options.duplicates, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if len(splitter) > 1:
			duplicates[splitter[0]] = splitter[1].split(';')

print("Assigning locations..")
add_data_from_data_and_duplicates_files(tree, countries_dict, options.duplicates, "country")

print("Collecting neighbours.. ")
entries = {}
with open(options.entry_nodes_file, "r") as enf:
	for line in enf:
		print(">entry_node " + line.strip())
		entry_node = tree.search_nodes(name=line.strip())[0]
		entry_strains = [n.name for n in entry_node.iter_leaves()]
		entry_countries = [countries_dict.get(s, "unknown") for s in entry_strains]
		entry_dates = [dates_dict.get(s, "") for s in entry_strains if countries_dict.get(s, "unknown") == "Russia"]
		for s in entry_strains:
			if s in duplicates:
				entry_countries.extend([countries_dict.get(d, "unknown") for d in duplicates[s]])
				if countries_dict.get(s, "unknown") == "Russia":
					entry_dates.extend([dates_dict.get(d, "") for d in duplicates[s]])
		entries[entry_node.name] = {"count":len(entry_countries),"countries_str":countries_collection_to_str(entry_countries), "max_date":max([d for d in entry_dates if not d == ""]), "min_date":min([d for d in entry_dates if not d == ""])}
		#entry_node.add_feature("strains_str", ",".join(entry_strains))
		#entry_node.add_feature("countries_str", countries_collection_to_str(entry_countries))

		farthest, maxdist = tree.get_farthest_node(line.strip())

		#  process closest foreign nodes
		output_foreign_nodes = []
		tnode = entry_node
		fmindist = maxdist

		while(not tnode.is_root() and get_distance_to_parent(child=entry_node, parent=tnode) <= fmindist):
			siss = tnode.get_sisters()
			for s in siss:
				curdist = get_distance_to_parent(child=entry_node, parent=tnode) + tnode.dist + s.dist
				fmindist, output_foreign_nodes = collectClosest(s, curdist, fmindist, output_foreign_nodes, leaf_not_fully_russian)
			tnode = tnode.up

		output_foreign_nodes.append(entry_node) # we need to keep it
		for n in output_foreign_nodes:
			n.add_features(keep=True)


print("Stripping tree..")
t2 = strip_tree(tree, options.output)
t2.write(format=1, outfile=options.output + "_test2")

print("Styling tree..")
add_data_from_data_and_duplicates_files(t2, countries_dict, options.duplicates, "country")
for e in entries:
	entry_node = t2.search_nodes(name=e)[0]
	#entry_node.add_feature("strains_str", entries[e]["strains_str"])
	entry_node.add_feature("countries_str", entries[e]["countries_str"])
	entry_node.add_feature("count", entries[e]["count"])
	entry_node.add_feature("min_date", entries[e]["min_date"])
	entry_node.add_feature("max_date", entries[e]["max_date"])
color_tree(t2)

ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = multiple_country_layout
ts.scale = 1
t2.ladderize()
#t2.render(options.output + ".png", tree_style=ts, dpi=300, w=870, units="mm")
t2.render(options.output + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
t2.render(options.output + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")


