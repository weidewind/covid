import pandas as pd
from ete3 import Tree
import re

merged = "/export/home/popova/workspace/covid/data/munched/gennady/definitely_all_nonrus_duplicates_of_nonrus_strains_3"
tree = "/export/home/popova/workspace/covid/data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted.rus_refluffed"
options_rpn_meta = "/export/home/popova/workspace/covid/data/russian/meta_all.xlsx"
options_meta = "/export/home/popova/workspace/covid/data/raw/metadata_2021-03-12_09-12.tsv"

print("Parsing metas..")
meta = pd.read_excel(options_rpn_meta)[['Внутренний номер']]
meta.columns = ['seq_id']
meta['country'] = 'Russia'
gismeta = pd.read_csv(options_meta, sep="\t")[['gisaid_epi_isl', 'country']]
gismeta.columns = ['seq_id', 'country']
print(gismeta.head())
meta = pd.concat([meta, gismeta])
meta_dict = dict(zip(meta['seq_id'], meta['country']))


print("Parsing tree..")
t = Tree(tree, format=1)
leaves = set([n.name for n in t.get_leaves()])
print("there are " + str(len(leaves)) + " in the tree\n")

print("Checking..")
with open(merged) as m:
	for line in m:
		splitter = line.strip().split("\t")
		tsplitter = re.split(';|\t',line.strip())
		if splitter[0] not in leaves:
			print("Wow, " + splitter[0] + " is not in the tree! that could only mean this is a duplicate of some russian strain, and it must be present in the rus_refluffed tree!")
		for strain in tsplitter:
			if meta_dict.get(strain) == "Russia":
				print("Oh my! " + strain + " is russian!")

