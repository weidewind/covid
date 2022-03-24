from meta import get_origin_dict
from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing origins..")
meta_dict = get_origin_dict()

print("Writing output..")
with open(options.output, "w") as out:
	for strain, origin in meta_dict.items():
		if not origin == "":
			out.write(strain + "\t" + origin + "\n")
	# 	out.write(strain + "\t" + origin + "\n")
	# for leaf in tree.iter_leaves():
	# 	if leaf.name not in meta_dict:
	# 		out.write(leaf.name + "\t" + "" + "\n")
		



