import pandas as pd
import re
import datetime
import optparse
import os.path
from ete3 import Tree

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-n', '--node', help='', type='str')
parser.add_option('-c', '--number', help='number of output trees', type='int')
parser.add_option('-s', '--step', help='number of nodes to go up', type='int')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)
node = tree.search_nodes(name=options.node)[0]

for i in range(1,options.number):
	node.write(outfile = options.output + "_" + str(i-1) + "_stepup" + str(options.step) + ".nwk", format=1)
	for j in range(1, options.step+1):
		node = node.up
