import pandas as pd
from ete3 import Tree
import re
from meta import get_country_dict
import optparse

parser = optparse.OptionParser()
parser.add_option('-m', '--merged', help='merged file for checking', type='str')
parser.add_option('-t', '--tree', help='', type='str')


options, args = parser.parse_args()

#merged = "/export/home/popova/workspace/covid/data/munched/gennady/definitely_all_nonrus_duplicates_of_nonrus_strains_3"
#tree = "/export/home/popova/workspace/covid/data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted.rus_refluffed"

meta_dict = get_country_dict()


print("Parsing tree..")
t = Tree(options.tree, format=1)
leaves = set([n.name for n in t.get_leaves()])
print("there are " + str(len(leaves)) + " in the tree\n")

print("Checking..")
with open(options.merged) as m:
	for line in m:
		splitter = line.strip().split("\t")
		tsplitter = re.split(';|\t',line.strip())
		if splitter[0] not in leaves:
			print("Wow, " + splitter[0] + " is not in the tree! that could only mean this is a duplicate of some russian strain, and it must be present in the rus_refluffed tree!")
		for strain in tsplitter:
			if meta_dict.get(strain) == "Russia":
				print("Oh my! " + strain + " is russian!")

