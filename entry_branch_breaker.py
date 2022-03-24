from ete3 import Tree, TreeNode
import time
import pandas as pd
import optparse
import sys
sys.setrecursionlimit(3500)


parser = optparse.OptionParser()
parser.add_option('--tree', help='tree file', type='str')
parser.add_option('--translin_file', help='$translin_file', type='str')


options, args = parser.parse_args()

def break_branches(tree, translin_file):
	broken_dict = {}
	with open(translin_file) as entries:
		for e in entries:
			entry = tree.search_nodes(name=e.strip())[0]
			dist = entry.dist
			if (dist > 0):
				parent = entry.up
				entry = entry.detach()
				middle = TreeNode(name = entry.name + "_middle")
				middle.add_child(child = entry, dist = dist/2)
				parent.add_child(child = middle, dist = dist/2)
				broken_dict[entry.name] = middle.name
			else:
				broken_dict[entry.name] = entry.name
	return(broken_dict)


def write_file_with_replacement(infile, outfile, replacement_dict):
	with open(infile, "r") as inf:
		with open(outfile, "w") as out:
			for line in inf:
				splitter = line.rstrip("\n").split("\t")
				if len(splitter) > 1 and splitter[1] in replacement_dict and not splitter[0] == splitter[1]: # change in .entries.exports, but in .out lines like "EPI_ISL01 EPI_ISL01" only change the first word
						splitter[1] = replacement_dict[splitter[1]]
				if len(splitter) == 3 and not splitter[2] == "":
					splitter[2] = replacement_dict[splitter[2]] # export nodes for secondary entries
				splitter[0] = replacement_dict[splitter[0]]
				newline = "\t".join(splitter)
				out.write(newline + "\n")


print("Parsing tree..")
tree = Tree(options.tree, format = 1)

print("Breaking entry branches.. ")
broken_dict_entries = break_branches(tree = tree, translin_file = options.translin_file + ".entries")
broken_dict_exports = break_branches(tree = tree, translin_file = options.translin_file + ".exports")
print(broken_dict_entries)
print(broken_dict_exports)
broken_dict = {**broken_dict_entries, **broken_dict_exports}  # merge two dicts (an entry cannot be an export, thus no conflicts expected)
print(broken_dict)

print("Writing broken tree..")
tree.write( format=1, outfile=options.tree + ".broken", format_root_node=True)

print("Writing broken entries and exports..")
write_file_with_replacement(options.translin_file + ".entries", options.translin_file + ".broken.entries", broken_dict)
write_file_with_replacement(options.translin_file + ".exports", options.translin_file + ".broken.exports", broken_dict)
write_file_with_replacement(options.translin_file, options.translin_file + ".broken", broken_dict)
write_file_with_replacement(options.translin_file + ".entries.exports", options.translin_file + ".broken.entries.exports", broken_dict)
