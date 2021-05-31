from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import time
import optparse
import matplotlib
import pandas as pd
import re



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s', '--states', help='tab-delimited file with states for all nodes', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

states = {}

# should this node be considered a leaf? 
def last_kept(node):
	if node.is_leaf():
		return(True)
	elif all(not ch.keep for ch in node.get_children()):
		return(True)
	else:
		return(False)

def strip_tree(t, outfile):
	# left, then right, then root
	for node in t.traverse("postorder"):
		if node.is_leaf():
			if states[node.name] == 1:
				node.add_features(keep=True)
			else:
				node.add_features(keep=False)
		if hasattr(node, 'keep') and node.keep == True:
			if not node.is_root():
				node.up.add_features(keep=True)
				for ch in node.up.get_children():
					ch.add_features(keep=True)
		else:
			node.add_features(keep=False)

	t.write(is_leaf_fn=last_kept, format=1, outfile=outfile, format_root_node=True)
	t2 = Tree( outfile, format=1)
	return (t2)


def add_countries_from_gismeta(t, meta_file):
	print("Reading gisaid metadata..")
	gismeta = pd.read_csv(meta_file, sep="\t")[['gisaid_epi_isl', 'country']]
	gismeta.columns = ['seq_id', 'region']
	print(gismeta.head())
	add_countries(t, gismeta)



def add_countries_from_countries_file(t, countries_file):
	countries = pd.read_csv(countries_file, sep="\t")
	countries.columns = ['seq_id', 'state', 'region']
	print(countries.head())
	add_countries(t, countries)


def add_countries(t, df):
		for node in t.traverse():
		if node.is_leaf() and re.match(r"^[a-zA-Z]",node.name): # after pruning some internal nodes became leaves, hence the additional check
			print(node.name)
			reg = df.loc[df['seq_id'] == node.name, "region"]
			print(reg)
			if reg.empty:
				node.add_features(country="unknown")
			else:
				node.add_features(country=reg.values[0])	


def color_tree(t):
	for node in t.traverse():
		print(node.name)
		nstyle = NodeStyle()
		nstyle["fgcolor"] = matplotlib.colors.to_hex([states[node.name], 0.0, 0.0])
		nstyle["size"] = 1
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		node.set_style(nstyle)


def country_layout(node):
	if node.is_leaf() and re.match(r"^[a-zA-Z]",node.name) and any(states[sis.name] > 0 for sis in node.up.get_children()):
		faces.add_face_to_node(AttrFace("country"), node, column=0)


print("Parsing tree..")
now = time.time()
tree = Tree(options.tree, format=1)
later = time.time()
diff = int(later - now)
print("Finished parsing in " + str(diff) + " seconds")


print("Parsing states..")
with open(options.states) as st:
	for line in st:
		[nname, _, state2_prob] = line.split()
		states[nname] = float(state2_prob)  # probability of being russian
		print("nname " + nname + " state 2 probability " + str(state2_prob))


t2 = strip_tree(tree, options.output)
t2.write(format=1, outfile=options.output + "test_t2")

add_countries_from_countries_file(t2, options.countries)

color_tree(t2)

ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = country_layout
ts.scale = 1
t2.ladderize()
#t2.render(options.output + ".png", tree_style=ts, dpi=300, w=870, units="mm")
t2.render(options.output + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
t2.render(options.output + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")
