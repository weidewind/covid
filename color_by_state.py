import optparse
import pandas as pd
import numpy as np
import re
from collections import Counter
from meta import get_region_dict
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace
import webcolors

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s', '--states', help='leaf_states.csv', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

states = pd.read_csv(options.states, sep="\s", header = None)
states_dict = dict(zip(states.iloc[:,0], states.iloc[:,2]))

print("Styling tree..")

def color_tree(t):
	for node in t.traverse():
		nstyle = NodeStyle()
		nstyle["size"] = 4
		nstyle["vt_line_width"] = 1
		nstyle["hz_line_width"] = 1
		nstyle["vt_line_color"] = "#aaaaaa"
		nstyle["hz_line_color"] = "#aaaaaa"
		if len(node.get_children()) == 1:
			node.add_face(TextFace(node.name, fgcolor="#000000"), column=0)
		prob = states_dict[node.name]
		print("#" + str(hex(int(prob*255))[2:].zfill(2)) + "0000")
		nstyle["fgcolor"] = "#" + str(hex(int(prob*255))[2:].zfill(2)) + "0000"
		# if prob == 0:
		# 	nstyle["fgcolor"] = "#000000"
		# if prob == 0.5:
		# 	nstyle["fgcolor"] = "#990000"
		# if prob == 1:
		# 	nstyle["fgcolor"] = "#ff0000"
		
		node.set_style(nstyle)

ts = TreeStyle()
ts.show_leaf_name = False

color_tree(tree)

tree.ladderize()
tree.render(options.output + ".svg", tree_style=ts, dpi=300, w=10000, units="mm")
tree.render(options.output + ".pdf", tree_style=ts, dpi=300, w=10000, units="mm")