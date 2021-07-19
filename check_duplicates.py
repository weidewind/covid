from ete3 import Tree
import pandas as pd
from meta import get_country_dict
import sys
import optparse

parser = optparse.OptionParser()
parser.add_option('-r', '--rtr_duplicates_file', help='_rus_to_rus_dups.tsv', type='str')
parser.add_option('-d', '--rtn_duplicates_file', help='_rus_to_nonrus_dups.tsv', type='str')
parser.add_option('-n', '--ntn_duplicates_file', help='_nonrus_to_nonrus_dups.tsv', type='str')
parser.add_option('-t', '--tree', help='_for_iqtree_3_-6.fasta.treefile', type='str')


options, args = parser.parse_args()


# check that all russian strains are russian, all non-russian strains are non-russian,
# and that all leader russian strains are present in the tree, and duplicates are not present in the tree

tree = Tree(options.tree, format = 1)
leaves_names = [n.name for n in tree.get_leaves()]

meta_dict = get_country_dict()


print("Parsing rus duplicates..")
with open(options.rtr_duplicates_file, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		if splitter[0] not in leaves_names:
			sys.stderr.write(splitter[0] + " main strain from rus to rus duplicates file " + options.rtr_duplicates_file + " is not in the tree!")
			raise
		for d in splitter[1].split(";"):
			if d in leaves_names:
				sys.stderr.write(d + " is a duplicate, but i found it in the tree!")
				#raise
		for s in splitter:
			dups = s.split(';')
			for d in dups:
				if meta_dict.get(d, "NA") != "Russia":
					sys.stderr.write(d + " from rus to rus duplicates file " + options.rtr_duplicates_file + " is not russian! " + meta_dict.get(d, "NA"))
					raise
					


print("Parsing rus to nonrus duplicates..")
with open(options.rtn_duplicates_file, "r") as dfile:
	for line in dfile:
		rus, nonrus = line.strip().split('\t')
		if rus not in leaves_names:
			sys.stderr.write(rus + " is not in the tree!")
			raise
		if meta_dict.get(rus, "NA") != "Russia":
				sys.stderr.write(rus + " main strain from rus to nonrus duplicates file " + options.rtn_duplicates_file + " is not russian! " + meta_dict.get(rd, "NA"))
				raise
		ndups = nonrus.split(';')
		for nd in ndups:
			if meta_dict.get(nd, "NA") == "Russia":
				sys.stderr.write(rd + " duplicate strain from rus to nonrus duplicates file " + options.rtn_duplicates_file + "  is  russian! " + meta_dict.get(nd, "NA"))
				raise
			if nd in leaves_names:
				sys.stderr.write(nd + " from rus to nonrus duplicates file " + options.rtn_duplicates_file + "  is a duplicate, but i found it in the tree!")
				#raise


print("Parsing nonrus duplicates..")
with open(options.ntn_duplicates_file, "r") as dfile:
	for line in dfile:
		splitter = line.strip().split('\t')
		for s in splitter:
			dups = s.split(';')
			for d in dups:
				if meta_dict.get(d, "NA") == "Russia":
					sys.stderr.write(d + " from nonrus to nonrus duplicates file " + options.ntn_duplicates_file + " is russian! " + meta_dict.get(d, "NA"))
					raise