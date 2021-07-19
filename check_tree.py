from ete3 import Tree
import pandas as pd


tree = "/export/home/popova/workspace/covid/data/munched/usher/nextstrain_and_17june_2/uncondensed-final-tree.nh"
meta = "/export/home/popova/workspace/covid/data/raw/17_06_2021_meta.csv"

tree = Tree(tree, format = 1)
leaves_names = [n.name for n in tree.get_leaves()]
print(leaves_names[0:5])

meta = pd.read_csv(meta, sep="\t")[['Внутренний номер', 'gisaid_id']]
rus_to_gis = dict(zip(meta['Внутренний номер'], meta['gisaid_id']))

print(dict(list(rus_to_gis.items())[0:5]))

for rus, gis in rus_to_gis.items():
	if rus in leaves_names and gis in leaves_names:
		print(rus + "," + gis)