from meta import get_country_dict
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import re



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()
print(options)

tree = Tree(options.tree, format = 1)
#print(tree.get_ascii(attributes=["name"], show_internal=True))
for n in tree.traverse("preorder"):
	if not n.is_leaf():
		pname = n.up.name if not n.is_root() else "sun itself"
		#print("dealing with " + n.name + ", whose parent name is " + pname + " and children names are " + ",".join([ch.name for ch in n.get_children()]))
		if len(n.get_children()) == 1:
			print(n.name + " has only one child!")		
			singleton = n.get_children()[0]
			singleton.dist = singleton.dist + n.dist
			if not n.name == "" and not singleton.is_leaf():
				singleton.name = n.name
			n.delete()
		#print("------")
#print(tree.get_ascii(attributes=["name"], show_internal=True))
tree.write(format=1, outfile=options.output)
