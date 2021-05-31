from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import time
from collections import Counter
import optparse
import matplotlib
import pandas as pd
import re


def add_data_from_data_and_duplicates_files(t, countries_dict, duplicates_file, feature_name):

	duplicates = {}
	with open(duplicates_file, "r") as dfile:
		for line in dfile:
			splitter = line.strip().split('\t')
			if len(splitter) > 1:
				duplicates[splitter[0]] = splitter[1].split(';')
	print(dict(list(duplicates.items())[0:5]))

	print("Assigning data to strains and duplicates..")
	for node in t.iter_leaves():
		print("Node " + node.name)
		if re.match(r"^[a-zA-Z]", node.name): # after pruning some internal nodes became leaves, hence the additional check
			node_countries = []
			reg = countries_dict[node.name]
			node_countries.append(reg)
			if node.name in duplicates:
				print(node.name + " has duplicates!")
				for dup in duplicates[node.name]:
					#print("found duplicate " + dup)
					if dup in countries_dict:
						node_countries.append(countries_dict[dup])
					else:
						node_countries.append("unknown")
			node.add_feature(feature_name, node_countries)
			#print("node "+node.name+" "+",".join(node_countries))