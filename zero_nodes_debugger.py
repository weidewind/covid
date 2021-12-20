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
parser.add_option('-z', '--zero', help='what length is considered zero', type='float', default = 0)
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()
print(options)

tree = Tree(options.tree, format = 1)
for n in tree.traverse():
	if n.name[:5] == "node_" and n.dist <= options.zero:
		print(n.name[:5] + " " + str(n.dist))
		n.delete()
tree.write(format=1, outfile=options.output, format_root_node=True)
