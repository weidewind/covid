#from Bio import Phylo
from ete3 import Tree
import time
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s', '--states', help='tab-delimited file with states for all nodes', type='str')
parser.add_option('-d', '--duplicates', help='file with duplicates (and cluster memebers)', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')
parser.add_option('-p', '--print_broader_subtree', default=-1, help='return newick subtrees for each entry. Make N steps up the tree before extracting a subtree', type='int')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

lineages = {} # lineages{entry_node_name} = @strains


def find_primary_lineages(node, array, entry_nname):
	array.append(node)
	exported = entry_nname is not None and states[node.name] == 0
	print("Processing node " + node.name + " with state 2 probability " + str(states[node.name]) + " and exproted = " + str(exported))
	#print("Is it exported? " + str(exported))
	#if len(lineages) > 0:
	#	print("Lineage len " + str(len(lineages)))
	if not exported:
		if entry_nname is not None:
			if node.is_leaf():
				lineages[entry_nname].append(node.name)
		else:
			if states[node.name] == 1:  # probability of being russian is 1
				entry_nname = node.name
				print("Lineage started at " + entry_nname)
				lineages[entry_nname] = []
				if node.is_leaf():
					lineages[entry_nname].append(node.name)
				#else:
					 #node.write(format=1, outfile=options.output + "_entry_" + node.name + ".newick")
	if not node.is_leaf() and not exported:
		for child in node.children:
			array = find_lineages(child, array, entry_nname)
	node = array.pop()
	return(array)


def find_lineages(node, array, entry_nname):
	array.append(node)
	exported = entry_nname is not None and states[node.name] == 0
	if exported:
		entry_nname = None
	print("Processing node " + node.name + " with state 2 probability " + str(states[node.name]) + " and exproted = " + str(exported))
	#print("Is it exported? " + str(exported))
	#if len(lineages) > 0:
	#	print("Lineage len " + str(len(lineages)))
	if not exported and entry_nname is not None:
		if node.is_leaf():
			lineages[entry_nname].append(node.name)
	if entry_nname is None:
		if states[node.name] == 1:  # probability of being russian is 1
			entry_nname = node.name
			print("Lineage started at " + entry_nname)
			lineages[entry_nname] = []
			if node.is_leaf():
				lineages[entry_nname].append(node.name)
				#else:
					 #node.write(format=1, outfile=options.output + "_entry_" + node.name + ".newick")
	if not node.is_leaf():
		for child in node.children:
			array = find_lineages(child, array, entry_nname)
	node = array.pop()
	return(array)


print("Parsing tree..")
now = time.time()
tree = Tree(options.tree, format = 1)
later = time.time()
diff = int(later - now)
print("Finished parsing in " + str(diff) + " seconds")


print("Parsing states..")
states = {}
with open(options.states) as st:
	for line in st:
		[nname, _, state2_prob] = line.split()
		states[nname] = float(state2_prob)  # probability of being russian
		print("nname " + nname + " state 2 probability " + str(state2_prob))


array = []
print("Looking for lineages.. ")
find_lineages(tree.get_tree_root(), array, None)

print("Parsing countries..")
countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
countries_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(countries_dict.items())[0:5]))

print("Adding duplicates..")
duplicates = {}
with open(options.duplicates) as dup:
	for line in dup:
		splitter =  line.strip().split("\t")
		if len(splitter) >1:
			for value in lineages.values():
				[centroid, strain_list_str] = splitter
				if centroid in value:
					rus_strain_list = [s for s in strain_list_str.split(";") if s in countries_dict and countries_dict[s] == "Russia"]
					value.extend(rus_strain_list)



linsizes = {}
for entry_nname in lineages.keys():
	size = len(lineages[entry_nname])
	if not size in linsizes:
		linsizes[size] = []
	linsizes[size].append(entry_nname)

with open(options.output + "stats", "w") as out:
	for size in sorted(linsizes.keys()):
		out.write("Found " + str(len(linsizes[size])) + " lineages of size " + str(size) + ": "  + ",".join(linsizes[size]) + "\n")

with open(options.output, "w") as out:	
	for entry_nname in lineages.keys():
		out.write(entry_nname + "\t")
		out.write(";".join(lineages[entry_nname]) + "\n")

with open(options.output + ".entries", "w") as out:
	for entry_nname in lineages.keys():
		out.write(entry_nname + "\n")

steps = options.print_broader_subtree
if steps > -1:
	for entry_nname in lineages.keys():
		entry_node = tree.search_nodes(name=entry_nname)[0]
		if not entry_node.is_leaf():
			counter = 0
			printing_node = entry_node
			while counter < steps:
				printing_node = printing_node.up
				counter += 1
			printing_node.write(format=1, outfile=options.output + "_entry_" + entry_node.name + "_" + str(steps) + "_steps_up.newick")

