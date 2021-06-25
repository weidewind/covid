import optparse
import pandas as pd
import re
from collections import Counter
from meta import get_region_dict
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
from tree_utils import add_data_from_data_and_duplicates_files


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-u', '--states', help='states file', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-a', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('-d', '--duplicates', help='states file', type='str')
parser.add_option('-l', '--lineages_list', help='lineages_of_interest.list', type='str')
parser.add_option('-p', '--pangostats_file', help='transmission_lineages.pangostats', type='str')
parser.add_option('-o', '--output', help='output', type='str')
parser.add_option('-s', '--steps', help='go n nodes up before printing', type='int')
parser.add_option('-g', '--pangolined', help='rus.pangolined.withduplicates', type='str')
parser.add_option('-e', '--entries_file', help='transmission_lineages.withduplicates.out', type='str')

options, args = parser.parse_args()





def color_tree(t):
	for node in t.traverse():
		nstyle = NodeStyle()
		nstyle["size"] = 1
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		if node.name in entries:
			nstyle["fgcolor"] = "#ee0000"
			nstyle["size"] = 4
		node.set_style(nstyle)


def countries_collection_to_str(c_collection):
	stats = dict(Counter(c_collection))
	sorted_stats = sorted(stats.items(), key=lambda x: x[1], reverse=True)
	string = ";".join([c[0] + " " + str(c[1]) for c in sorted_stats])
	return (string)

entries = []
ilinlist = pd.read_csv(options.lineages_list)
ilins = ilinlist["lineage"]
print(ilins)

pangolist = pd.read_csv(options.pangostats_file, sep=',')
with open(options.output + ".entries", "w") as out:
	for index, row in pangolist.iterrows():
		for lin in row["lineage_stat"].split(";"):
			lin = lin.split(":")[0]
			for ilin in ilins:
				if lin == ilin or (len(lin) > len(ilin) and lin[:len(ilin)+1] == ilin + "."):
					print(" lin, ilin: " + lin + ", " + ilin)
					entries.append(row["entry"])
					out.write(row["entry"] + "\t" + lin + "\n")

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing countries..")
countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
countries_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(countries_dict.items())[0:5]))

print("Parsing regions..")
regions_dict = get_region_dict()
print(dict(list(regions_dict.items())[0:5]))

print("Parsing dates..")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))


print("Parsing states..")
states = {}
with open(options.states) as st:
	for line in st:
		[nname, _, state2_prob] = line.split()
		states[nname] = float(state2_prob)  # probability of being russian
		#print("nname " + nname + " state 2 probability " + str(state2_prob))

print("Parsing entries..")
entry_contents = {}
with open(options.entries_file) as ef:
	for line in ef:
		e, content = line.strip().split("\t")
		entry_contents[e] = content.split(";")

print("Parsing duplicates..")
duplicates = {}
with open(options.duplicates, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if len(splitter) > 1:
			duplicates[splitter[0]] = splitter[1].split(';')

print("Parsing pangolineages..")
pango = pd.read_csv(options.pangolined, sep=",")
pango_dict = dict(zip(pango['taxon'], pango['lineage']))
print(dict(list(pango_dict.items())[0:5]))

pango_color_dict = {"B.1.1.7":"#1b9e77", "B.1.351":"#d95f02", "B.1.617":"#7570b3", "AT.1":"#e7298a", "B.1.1.451":"#66a61e", "B.1.1.317":"#e6ab02"}



print("Drawing trees..")

for entry in entries:
	print("Entry " + entry)
	entry_node = tree.search_nodes(name=entry)[0]

	step = 0
	tnode = entry_node
	parents = []
	while step < options.steps:
		tnode = tnode.up
		parents.append(tnode.name)
		step += 1

	def is_foreign_child(node):
		if node.is_leaf():
			return(True)
		if states[node.name] == 0 and node.name not in parents:
			return(True)
		else:
			return(False)

	filename = options.output + "_entry_" + entry + ".newick"
	tnode.write(is_leaf_fn=is_foreign_child, format=1, outfile=filename, format_root_node=True)
	#tnode.write(format=1, outfile=filename)
	t = Tree(filename, format=1)

	print("Assigning locations..")
	add_data_from_data_and_duplicates_files(t, countries_dict, options.duplicates, "country")

	entry_node = t.search_nodes(name=entry)[0]

	#entry_strains = [n.name for n in entry_node.iter_leaves()]
	entry_strains = entry_contents[entry_node.name]
	#entry_countries = [countries_dict.get(s, "unknown") for s in entry_strains]
	entry_dates = [dates_dict.get(s, "") for s in entry_strains if countries_dict.get(s, "unknown") == "Russia"]
	for s in entry_strains:
		if s in duplicates:
			#entry_countries.extend([countries_dict.get(d, "unknown") for d in duplicates[s]])
			if countries_dict.get(s, "unknown") == "Russia":
				entry_dates.extend([dates_dict.get(d, "") for d in duplicates[s]])


	print("Styling tree..")

	entry_node.add_feature("count", len(entry_dates))
	#entry_node.add_feature("countries_str", countries_collection_to_str(entry_countries))
	entry_node.add_feature("min_date", min([d for d in entry_dates if not d == "" and not d == "unknown"]))
	entry_node.add_feature("max_date", max([d for d in entry_dates if not d == "" and not d == "unknown"]))
	color_tree(t)

	def multiple_country_layout(node):
		if node.is_leaf() and  re.match(r"^[a-zA-Z]", node.name) and countries_dict[node.name] == "Russia" and not node.name == entry_node.name:
			#faces.add_face_to_node(faces.TextFace(node.date[:11], fgcolor="#ee0000"), node, column=0)
			print(node.name)
			linsplitter = pango_dict[node.name].split(".")
			linup = ".".join(linsplitter[:-1])
			pango_color = pango_color_dict.get(pango_dict[node.name], pango_color_dict.get(linup, "#999999"))
			faces.add_face_to_node(faces.TextFace(node.name  + " " + dates_dict.get(node.name, "unknown")[:11] + " " + pango_dict[node.name] + " " + regions_dict.get(node.name, "unknown"), fgcolor=pango_color), node, column=0)
		elif node.is_leaf()  and re.match(r"^[a-zA-Z]", node.name):
			faces.add_face_to_node(faces.TextFace(node.name + " " + str(countries_collection_to_str(node.country))), node, column=0)
		#elif not node.name == entry_node.name:
		#	faces.add_face_to_node(faces.TextFace(node.name + " state:" + str(states[node.name])), node, column=0)
		if node.name == entry_node.name:
			faces.add_face_to_node(faces.TextFace(node.name, fgcolor="#ee0000"), node, column=0)
			if node.count > 1:
				faces.add_face_to_node(faces.TextFace("earliest " + node.min_date[:11], fgcolor="#ee0000"), node, column=0)
				faces.add_face_to_node(faces.TextFace("latest " + node.max_date[:11], fgcolor="#ee0000"), node, column=0)

	ts = TreeStyle()
	ts.show_leaf_name = False
	ts.layout_fn = multiple_country_layout
	ts.scale = 100000
	t.ladderize()
	t.render(options.output + "_entry_" + entry + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
	t.render(options.output + "_entry_" + entry + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")






