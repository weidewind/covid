from meta import get_country_dict, get_district_dict
from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Attributing to districts..")
meta_dict = get_district_dict()

print("Writing output..")
with open(options.output, "w") as out:
	for strain, region in meta_dict.items():
		out.write(strain + "\t" + region + "\n")
	for leaf in tree.iter_leaves():
		if leaf.name not in meta_dict:
			out.write(leaf.name + "\t" + "unknown" + "\n")
		



