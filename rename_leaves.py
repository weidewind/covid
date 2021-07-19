from ete3 import Tree
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-m', '--meta',type='str')
parser.add_option('-l', '--long_ids',help='file for ids conversion', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()

tree = Tree(options.tree, format=1)
meta = pd.read_csv(options.meta,sep='\t')[["Внутренний номер", "gisaid_id"]]
meta.columns = ["id", "gisaid_id"]
gisaid_dict = dict(zip(meta["gisaid_id"], meta["id"]))

long_ids = pd.read_csv(options.long_ids,sep='\t')
long_ids.columns = ["id", "long_id"]
long_id_dict = dict(zip(long_ids["long_id"], long_ids["id"]))

for leaf in tree.get_leaves():
	newname = leaf.name
	if newname in long_id_dict:
		newname = long_id_dict[newname]
	if newname in gisaid_dict:
		newname = gisaid_dict[newname]
	leaf.name = newname

print("Writing new tree to " + options.output)
tree.write(outfile=options.output, format=1)
