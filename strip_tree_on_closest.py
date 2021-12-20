from meta import get_country_dict
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import re



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-m', '--melted_foreign', help='covid/output/variant_analysis/../melted_foreign from get_melted_foreign_strains_for_entries.R used in select_rus_from_mmsa.sh', type='str')
parser.add_option('-d', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('-e', '--entry_nodes_file', help='file with entry nodes (one per line, output from find_transmission_lineages.py)', type='str')
#parser.add_option('-c', '--entry_mut_stats', help='file with counts of mutated and non-mutated strains for each entry', type='str')
#parser.add_option('-n', '--best_entries_number', help='number of best entries to highlight', type='int')
parser.add_option('-s', '--scale', help='pixels to tree length unit', type='int', default = 1)
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()
print(options)


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


def notRussian(leaf):
	if country_dict.get(leaf.name, "unknown") == "Russia":
		return(False)
	else:
		return(True)


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


def strip_all_but_entries(t, entries,outfile):
	for node in t.traverse("postorder"):
		if node.name in entries:
			node.add_features(keep=True)
		if hasattr(node, 'keep') and node.keep:
			if not node.is_root():
				node.up.add_features(keep=True)
		else:
			node.add_features(keep=False)
	t.write(is_leaf_fn=last_kept, format=1, outfile=outfile, format_root_node=True)
	t2 = Tree(outfile, format=1)
	return (t2)



# def strip_all_but_entries(t, entries):
# 	for node in t.traverse():
# 		if node.name not in entries and not node.is_root():
# 			node.delete()
# 			print("deleted " + node.name)
# 	return(t)

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


def mut_count_to_str(node):
	r = ""
	if node.mut_count > 0:
		r = "mut: " + str(node.mut_count)
	if node.nonmut_count > 0:
		r = r + " nonmut: " + str(node.nonmut_count)
	return (r)

def clade_mut_count_to_str(node):
	r = ""
	if node.clade_mut_count > 0:
		r = "clade_rus_mut: " + str(node.clade_mut_count)
	if node.nonmut_count > 0:
		r = r + " clade_rus_nonmut: " + str(node.clade_nonmut_count)
	return (r)	


def multiple_country_layout(node):
	if node.is_leaf() and re.match(r"^[a-zA-Z]", node.name) and node.name not in entries:
		if node.name in foreign_mut_dict:
			color =  "#000077" if foreign_mut_dict[node.name] == 1 else "#770000"
			faces.add_face_to_node(faces.TextFace(country_dict.get(node.name, "unknown")  + " " + dates_dict.get(node.name, "unknown"), fgcolor=color), node, column=0)
		else:
			faces.add_face_to_node(faces.TextFace(country_dict.get(node.name, "unknown")), node, column=0)
	if node.name in entries:
		color =  "#0000ee" if node.mut_count > node.nonmut_count else "#ee0000"
		faces.add_face_to_node(faces.TextFace(node.name + " " + clade_mut_count_to_str(node) + " " + mut_count_to_str(node) + " " + str(node.countries_str) , fgcolor=color), node,  column=0)
		if node.mut_count + node.nonmut_count  > 1:
			faces.add_face_to_node(faces.TextFace("earliest " + node.min_date[:11], fgcolor=color), node,  column=0)
			faces.add_face_to_node(faces.TextFace("latest " + node.max_date[:11], fgcolor=color), node,  column=0)
		else:
			faces.add_face_to_node(faces.TextFace(node.min_date[:11], fgcolor=color), node,  column=0)
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


def entries_layout(node):
	if node.name in entries:
		color =  "#0000ee" if node.mut_count > node.nonmut_count else "#ee0000"
		faces.add_face_to_node(faces.TextFace(node.name + " " + clade_mut_count_to_str(node) + " " + mut_count_to_str(node), fgcolor=color), node,  column=0)
		if node.mut_count + node.nonmut_count  > 1:
			faces.add_face_to_node(faces.TextFace("earliest " + node.min_date[:11], fgcolor=color), node,  column=0)
			faces.add_face_to_node(faces.TextFace("latest " + node.max_date[:11], fgcolor=color), node,  column=0)
		else:
			faces.add_face_to_node(faces.TextFace(node.min_date[:11], fgcolor=color), node,  column=0)

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

# best_entries = []
# if options.entry_mut_stats:
# 	estats = pd.read_csv(options.entry_mut_stats)
# 	best_entries = estats[0:options.best_entries_number]["entry"]
# 	print("best entries: ")
# 	print(best_entries)

print("Collecting entries info.. ")
entries = {}
variants = pd.read_csv(options.entry_nodes_file, sep = ",",  dtype={'entry': str})
rus_mut_dict = dict(zip(variants['taxon'], variants['is_double_mut']))
unique_entries = list(set(variants["entry"]))

print("Parsing dates and countries..")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))

country_dict = get_country_dict()
print(dict(list(country_dict.items())[0:5]))

print("Collecting neighbours.. ")
f = pd.read_csv(options.melted_foreign, sep = ",")
print(f[0:5])
foreign_mut_dict = dict(zip(f['taxon'], f['is_double_mut']))
melted_foreign_strains = f["taxon"]




for e in unique_entries:
	print(">entry_node " + e)
	entry_node = tree.search_nodes(name=e)[0]
	entry_node.add_features(keep=True)
	entry_strains = [n.name for n in entry_node.iter_leaves()]
	clade_mut_count = len([s for s in entry_strains if rus_mut_dict.get(s, "unknown") == 1])
	clade_nonmut_count = len([s for s in entry_strains if rus_mut_dict.get(s, "unknown") == 0])
	entry_countries = [country_dict.get(s, "unknown") for s in entry_strains]
	entry_dates = [dates_dict.get(s, "") for s in entry_strains if country_dict.get(s, "unknown") == "Russia"]
	mut_count = variants[(variants["entry"] == e) & (variants["is_double_mut"] == 1)]["taxon"].count()
	nonmut_count = variants[(variants["entry"] == e) & (variants["is_double_mut"] == 0)]["taxon"].count()
	entries[entry_node.name] = {"count":len(entry_countries),"clade_mut_count":clade_mut_count, "clade_nonmut_count":clade_nonmut_count,"mut_count":mut_count,"nonmut_count":nonmut_count,"countries_str":countries_collection_to_str(entry_countries), "max_date":max([d for d in entry_dates if not d == ""]), "min_date":min([d for d in entry_dates if not d == ""])}

	farthest, maxdist = tree.get_farthest_node(e)

	#  process closest foreign nodes
	output_foreign_nodes = []
	tnode = entry_node
	fmindist = maxdist

	while(not tnode.is_root() and get_distance_to_parent(child=entry_node, parent=tnode) <= fmindist):
		siss = tnode.get_sisters()
		for s in siss:
			curdist = get_distance_to_parent(child=entry_node, parent=tnode) + tnode.dist + s.dist
			fmindist, output_foreign_nodes = collectClosest(s, curdist, fmindist, output_foreign_nodes, notRussian)
		tnode = tnode.up
	print("found " + str(len(output_foreign_nodes)) + " closest foreign nodes")
	for n in output_foreign_nodes:
		if n.name in melted_foreign_strains:
			n.add_features(keep=True)


print("Stripping tree..")
t2 = strip_tree(tree, options.output)
t2.write(format=1, outfile=options.output + ".nw")

print("Styling tree..")
for e in entries:
	entry_node = t2.search_nodes(name=e)[0]
	#entry_node.add_feature("strains_str", entries[e]["strains_str"])
	entry_node.add_feature("countries_str", entries[e]["countries_str"])
	#entry_node.add_feature("count", entries[e]["count"])
	entry_node.add_feature("mut_count", entries[e]["mut_count"])
	entry_node.add_feature("nonmut_count", entries[e]["nonmut_count"])
	entry_node.add_feature("clade_mut_count", entries[e]["clade_mut_count"])
	entry_node.add_feature("clade_nonmut_count", entries[e]["clade_nonmut_count"])
	entry_node.add_feature("min_date", entries[e]["min_date"])
	entry_node.add_feature("max_date", entries[e]["max_date"])
color_tree(t2)

ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = multiple_country_layout
ts.scale = options.scale
t2.ladderize()
#t2.render(options.output + ".png", tree_style=ts, dpi=300, w=870, units="mm")
t2.render(options.output + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
t2.render(options.output + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")

print("Stripping all but entries..")
tree = Tree(options.tree, format=1)
t3 = strip_all_but_entries(tree, unique_entries,  options.output +"_only_entries.nw")
print("Styling tree..")
for e in entries:
	print("entry " + e)
	entry_node = t3.search_nodes(name=e)[0]
	#entry_node.add_feature("strains_str", entries[e]["strains_str"])
	entry_node.add_feature("countries_str", entries[e]["countries_str"])
	#entry_node.add_feature("count", entries[e]["count"])
	entry_node.add_feature("mut_count", entries[e]["mut_count"])
	entry_node.add_feature("nonmut_count", entries[e]["nonmut_count"])
	entry_node.add_feature("clade_mut_count", entries[e]["clade_mut_count"])
	entry_node.add_feature("clade_nonmut_count", entries[e]["clade_nonmut_count"])
	entry_node.add_feature("min_date", entries[e]["min_date"])
	entry_node.add_feature("max_date", entries[e]["max_date"])
color_tree(t3)
ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = entries_layout
ts.scale = options.scale
t3.ladderize()
#t2.render(options.output + ".png", tree_style=ts, dpi=300, w=870, units="mm")
t3.render(options.output + "_only_entries.svg", tree_style=ts, dpi=300, w=1000, units="mm")
t3.render(options.output + "_only_entries.pdf", tree_style=ts, dpi=300, w=1000, units="mm")
