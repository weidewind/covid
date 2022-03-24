from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s', '--countries', help='leaf_states.csv', type='str')
parser.add_option('-d', '--pangolined', help='file with pangolin annotation', type='str')
parser.add_option('-o', '--output', help='output tree', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)
leaves = tree.get_leaves()

if options.countries:
	print("Parsing countries..")
	countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
	countries_dict = dict(zip(countries['seq_id'], countries['region']))
	print(dict(list(countries_dict.items())[0:5]))

pangodict = {}
pangofiles = options.pangolined.split(",")
for pfile in [f for f in pangofiles if f != ""]:
	with open(pfile, "r") as pango:
		for l in pango:
			splitter = l.strip().split(",")
			pangodict[splitter[0]] = ",".join(splitter[1:])

with open(options.output, "w") as out:
	for l in leaves:
		if l.name not in pangodict and (not options.countries or countries_dict.get(l.name, "") == "Russia"):
			out.write(l.name + "\n")