from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--pangolined', help='file with pangolin annotation', type='str')
parser.add_option('-o', '--output', help='output tree', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)
leaves = tree.get_leaves()

pangodict = {}
with open(options.pangolined, "r") as pango:
	for l in pango:
		splitter = l.strip().split(",")
		pangodict[splitter[0]] = ",".join(splitter[1:])

with open(options.output, "w") as out:
	for l in leaves:
		if l.name not in pangodict:
			out.write(l.name + "\n")