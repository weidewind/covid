import pandas as pd
import re
import datetime
import optparse
import os.path
from ete3 import Tree

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-p', '--pangolined', help='', type='str')
parser.add_option('-s', '--states', help='c for circular or r for rectangular', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)
leaves = tree.get_leaves()

print("Parsing states..")
print(options.states)
countries = pd.read_csv(options.states, sep="\t", names=['seq_id', 'state', 'region'])
countries_dict = dict(zip(countries['seq_id'], countries['state']))
print(dict(list(countries_dict.items())[0:5]))

print("Parsing pangolineages..")
pango_dict = {}
header = ""
with open(options.pangolined, "r") as pango:
	header = pango.readline()
	print("header " + header)
	for l in pango:
		splitter = l.strip().split(",")
		pango_dict[splitter[0]] = l

check_dupls = {}
print("Writing output..")
with open(options.output, "w") as out:
	out.write(header)
	for l in leaves:
		if countries_dict.get(l.name, "") == 2 and l.name not in check_dupls:
			check_dupls[l.name] = 1
			out.write(pango_dict[l.name])

