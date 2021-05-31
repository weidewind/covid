from meta import get_country_dict
from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-m', '--meta', help='csv meta file form gisaid', type='str')
parser.add_option('-r', '--rpn_meta', help='xlsx rpn meta file', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

meta_dict = get_country_dict()

print("Writing output..")
with open(options.output, "w") as out:
	for strain, country in meta_dict.items():
		state = 2 if country == "Russia" else 1
		out.write(strain + "\t" + str(state) + "\t" + country + "\n")
	for leaf in tree.iter_leaves():
		if leaf.name not in meta_dict:
			out.write(leaf.name + "\t" + str(1) + "\t" + "unknown" + "\n")
		
