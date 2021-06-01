from ete3 import Tree
import pandas as pd
from meta import get_country_dict

# check that all russian strains are russian, all non-russian strains are non-russian,
# and that all  leader russian strains are present in the tree

tree_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_for_iqtree_3_-6.fasta.treefile"
duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_rus_to_rus_dups.tsv"
rtn_duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_rus_to_nonrus_dups.tsv"
ntn_duplicates_file = "/export/home/popova/workspace/covid/data/munched/gennady/may31/21_05_2021_nonrus_to_nonrus_dups.tsv"

tree = Tree(tree_file, format = 1)
leaves_names = [n.name for n in tree.get_leaves()]

meta_dict = get_country_dict()


print("Parsing rus duplicates..")
with open(duplicates_file, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if splitter[0] not in leaves_names:
			print(splitter[0] + " is not in the tree!")
		for d in splitter[1].split(";"):
			if d in leaves_names:
				print(d + " is a duplicate, but i found it in the tree!")
		for s in splitter:
			dups = s.split(';')
			for d in dups:
				if meta_dict.get(d, "NA") != "Russia":
					print(d + " is not russian! " + meta_dict.get(d, "NA"))


print("Parsing rus to nonrus duplicates..")
with open(rtn_duplicates_file, "r") as dfile:
	for line in dfile:
		rus, nonrus = line.strip().split('\t')
		if rus not in leaves_names:
			print(rus + " is not in the tree!")
		if meta_dict.get(rus, "NA") != "Russia":
				print(rus + " is not russian! " + meta_dict.get(rd, "NA"))
		ndups = nonrus.split(';')
		for nd in ndups:
			if meta_dict.get(nd, "NA") == "Russia":
				print(rd + " is  russian! " + meta_dict.get(nd, "NA"))
			if nd in leaves_names:
				print(nd + " is a duplicate, but i found it in the tree!")


print("Parsing nonrus duplicates..")
with open(ntn_duplicates_file, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		for s in splitter:
			dups = s.split(';')
			for d in dups:
				if meta_dict.get(d, "NA") == "Russia":
					print(d + " is russian! " + meta_dict.get(d, "NA"))