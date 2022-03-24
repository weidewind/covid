import optparse
import pandas as pd
import numpy as np
import re
from collections import Counter
from meta import get_region_dict
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace
from tree_utils import add_data_from_data_and_duplicates_files
from datetime import datetime
import webcolors
from meta import get_date_dict
import random


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--dates', help='leaf_dates.csv', type='str')
parser.add_option('-s', '--states', help='sankoff probs', type='str')
parser.add_option('-c', '--countries', help='leaf_states.csv', type='str')
parser.add_option('--districts', help='leaf_districts.csv', type='str')
parser.add_option('-p', '--pangolined', help='', type='str')
parser.add_option('-l', '--lineage', help='lineage of interest', type='str')
parser.add_option('-n', '--node', help='node of interest', type='str')
parser.add_option('-m', '--mode', help='c for circular or r for rectangular', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

def style_tree(t):
	for node in t.traverse():
		nstyle = NodeStyle()
		nstyle["fgcolor"] = shade_of_grey(state_dict[node.name])
		nstyle["vt_line_width"] = 2
		nstyle["hz_line_width"] = 2
		nstyle["vt_line_color"] = shade_of_grey(state_dict[node.name])
		nstyle["hz_line_color"] = shade_of_grey(state_dict[node.name])
		if node.name == options.node:
			node.add_face(TextFace(node.name, fgcolor="#000000"), column=0) #fsize = size
		if node.is_leaf() and pango_dict[node.name][0:len(options.lineage)] == options.lineage:
			nstyle["size"] = 5
			nstyle["fgcolor"] = "limegreen"
		if node.is_leaf() and not dates_dict[node.name] == "unknown":
			date = datetime.strptime(dates_dict[node.name], "%Y-%m-%d")
			if int(dates_dict[node.name].split("-")[2])%10 == 0 or date < datetime.strptime("2021-03-10", "%Y-%m-%d") or node.up.name in [options.node, "265512", "265511"] or (not node.up.is_root() and node.up.up.name in [options.node, "265512", "265511"]):
				text = dates_dict[node.name] +  " " + node.name
				if not pango_dict[node.name] == "B.1.1.523":
						text = text + " " + pango_dict[node.name]
				if state_dict[node.name] == 1:
						text = text + " " + district_dict[node.name]
				else:
						text = text + " " + country_dict[node.name]

				if date < datetime.strptime("2021-03-10", "%Y-%m-%d"):
					fgc = "red"
				else:
					fgc = "black"
				node.add_face(TextFace(text,  fgcolor=fgc), column=0) #fsize = size
		node.set_style(nstyle)


def shade_of_grey(prob):
	if prob < 0.1:
		prob = 0.1
	t = str(hex(int(256*(1-prob)))[2:].zfill(2))
	return("#"+t+t+t)

print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Parsing pangolineages..")
pango = pd.read_csv(options.pangolined, sep=",")
pango_dict = dict(zip(pango['taxon'], pango['lineage']))
print(dict(list(pango_dict.items())[0:5]))

print("Parsing countries..")
countries = pd.read_csv(options.countries, sep="\t", names=['seq_id', 'state', 'region'])
country_dict = dict(zip(countries['seq_id'], countries['region']))
print(dict(list(country_dict.items())[0:5]))

print("Parsing districts..")
districts = pd.read_csv(options.districts, sep="\t", names=['seq_id', 'region'])
district_dict = dict(zip(districts['seq_id'], districts['region']))
print(dict(list(district_dict.items())[0:5]))

print("Parsing states..")
states = pd.read_csv(options.states, sep="\s", names=['seq_id', 'foreign_prob', 'rus_prob'])
state_dict = dict(zip(states['seq_id'], states['rus_prob']))
print(dict(list(state_dict.items())[0:5]))

print("Parsing dates..")
dates_dict = get_date_dict(options.dates)
print(dict(list(dates_dict.items())[0:5]))


print("Styling tree..")

ts = TreeStyle()
ts.mode = options.mode
#ts.root_opening_factor = 1
ts.show_leaf_name = False
if options.mode == "r":
	ts.scale = 10

ts.allow_face_overlap = True
style_tree(tree)


tree.ladderize()
tree.render(options.output + ".svg", tree_style=ts, dpi=72, w=10000, units="mm")
tree.render(options.output + ".pdf", tree_style=ts, dpi=72, w=10000, units="mm")