from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--duplicates', help='file with duplicates', type='str')
parser.add_option('-o', '--output', help='output tree', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing duplicates..")
duplicates = {}
with open(options.duplicates, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if len(splitter) > 1:
			duplicates[splitter[0]] = splitter[1].split(';')

print("Adding duplicate leaves..")
for leaf in tree.get_leaves():
	if leaf.name in duplicates:
		leafn = leaf.name
		leaf.name = "node_" + leaf.name
		leaf.add_child(name=leafn, dist=0)
		for dup in duplicates[leafn]:
			leaf.add_child(name=dup, dist=0)

print("Writing new tree..")
tree.write(format=1, outfile=options.output)