from meta import get_date_dict
import datetime
from ete3 import Tree
import pandas as pd
import optparse
import re


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Making dictionary..")
meta_dict = get_date_dict()

print("Writing output..")
with open(options.output, "w") as out:
	for strain, date in meta_dict.items():
		out.write(strain + "\t" + date + "\n")
	for leaf in tree.iter_leaves():
		if meta_dict[leaf.name] =="unknown":
			print("Unknown date for tree leaf " + leaf.name)
		if leaf.name not in meta_dict:
			print("Unknown date for tree leaf " + leaf.name)
			out.write(leaf.name + "\t" + "unknown" + "\n")
