from tree_utils import add_data_from_data_and_duplicates_files
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import datetime
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import numpy as np
import re
import sys
import random

sys.setrecursionlimit(3500)

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-e', '--nodes_file', help='file with nodes', type='str', default=None)
parser.add_option('-o', '--output', help='output', type='str')
parser.add_option('-m', '--maxdist', help='maximum distance for closest nodes', type='float', default = None)
parser.add_option('-r', '--random', help='pick n random nodes', type='int', default = None)
parser.add_option('--shift', dest='shift', action='store_true', default=False)


options, args = parser.parse_args()


def collectClosestIngroup(entry_node, maxdist=None):
	output_nodes = []
	if maxdist is None:
		farthest, maxdist = tree.get_farthest_node()
	mindist = maxdist

	curdist = entry_node.dist/2 if options.shift else 0
	mindist, output_nodes = collectClosest(entry_node, curdist, mindist, output_nodes)

	return output_nodes, mindist


def collectClosestOutgroup(entry_node, maxdist=None):
	output_nodes = []
	if maxdist is None:
		farthest, maxdist = tree.get_farthest_node()
	mindist = maxdist
	tnode = entry_node

	while(not tnode.is_root() and get_distance_to_parent(child=entry_node, parent=tnode) <= mindist):
		siss = tnode.get_sisters()
		for s in siss:
			#print("sis " + s.name)
			shift = entry_node.dist/2 if options.shift else 0
			curdist = get_distance_to_parent(child=entry_node, parent=tnode) - shift + tnode.dist + s.dist
			mindist, output_nodes = collectClosest(s, curdist, mindist, output_nodes)
		tnode = tnode.up

	return output_nodes, mindist


def collectClosest(node, curdist, mindist, output_nodes):  # curdepth - current distance to the leaf in question; mindist - best distance so far
	#print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_nodes", ",".join([n.name for n in output_nodes])]))
	if node.is_leaf() and curdist <= mindist:
		if curdist < mindist:
			output_nodes.clear()
		output_nodes.append(node)
		mindist = curdist
	elif curdist <= mindist:
		for ch in node.children:
			chdist = curdist + ch.dist
			if chdist <= mindist: ## add or output_nodes does not contain dates
				(mindist, output_nodes) = collectClosest(ch, chdist, mindist, output_nodes)
		#print("After traversing all children of node " + node.name + ":")
		#print(" ".join(["node", node.name, "curdist", str(curdist), "mindist", str(mindist), "output_nodes", ",".join([n.name for n in output_nodes])]))
	return mindist, output_nodes


def get_distance_to_parent(child, parent):
	temp = child
	dist = 0
	while not temp.name == parent.name:
		dist += temp.dist
		temp = temp.up
	return(dist)


def process_node(entry_node, out):
	innodes, imindist = collectClosestIngroup(entry_node,  maxdist=options.maxdist)
	outnodes, omindist = collectClosestOutgroup(entry_node, maxdist=options.maxdist)
	
	out.write("\t".join([",".join([n.name for n in innodes]), str(imindist)])+ "\t")
	out.write("\t".join([",".join([n.name for n in outnodes]), str(omindist)])+ "\n")



print("Parsing tree..")
tree = Tree(options.tree, format=1)

print("Processing nodes..")
with open(options.output + ".dates_stats", "w") as out:
	out.write("\t".join(["entry", "closest_ingroup_nodes","tree_distance_to_closest_rus_strain", "closest_outgroup_nodes","tree_distance_to_closest_foreign_strain"]) + "\n")

	if not options.nodes_file is None:
		with open(options.nodes_file, "r") as nf:
			for line in nf:
				print("Collecting closest nodes..")
				print(">node " + line.strip())
				out.write(line.strip() + "\t")
				entry_node = tree.search_nodes(name=line.strip())[0]
				entry_dist_dict = {}
				process_node(entry_node, out)
	elif not options.random is None:
		nodes = [n for n in tree.traverse("postorder")]
		entry_nodes = random.sample(nodes, options.random)
		for entry_node in entry_nodes:
			if not (len(entry_node.name)>len("_middle") and entry_node.name[-7:] == "_middle"):
				print("Collecting closest nodes..")
				print(">node " + entry_node.name)
				out.write(entry_node.name + "\t")
				entry_dist_dict = {}
				process_node(entry_node, out)



			




