from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--duplicates', help='file with rus no nonrus duplicates', type='str')
parser.add_option('-o', '--output', help='output tree', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing duplicates..")
duplicates = {}
with open(options.duplicates, "r") as dfile:
	with open(options.output + ".duplicates", "w") as out:
		for line in dfile:
			splitter = line.strip().split('\t')
			if len(splitter) > 1:
				nonrus = splitter[1].split(';')
				duplicates[nonrus[0]] = nonrus[1:]
				print("removing duplicates of " + nonrus[0])
				node = tree.search_nodes(name=nonrus[0])[0]
				for n in node.up.get_children():
					if n.name in nonrus[1:]:
						n.detach()
				out.write(nonrus[0] + "\t" + ";".join(nonrus[1:]) + "\n")


print("Writing new tree..")
tree.write(format=1, outfile=options.output + ".newick")