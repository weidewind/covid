from ete3 import Tree
import time
import datetime
import pandas as pd
import optparse


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-d', '--dates', help='tab-delimited file with states for all nodes', type='str')
parser.add_option('-l', '--from_date', help='yyyy-mm-dd', type='str')
parser.add_option('-u', '--to_date', help='yyyy-mm-dd', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


def excess_leaf(node, from_date, to_date):
	if dates_dict[node.name] == "unknown":
		return True
	nodedate = datetime.datetime.strptime(dates_dict[node.name], '%Y-%m-%d') 
	if nodedate >= datetime.datetime.strptime(to_date, '%Y-%m-%d') or nodedate < datetime.datetime.strptime(from_date, '%Y-%m-%d'):
		return True
	else:
		return False


print("Parsing tree..")
tree = Tree(options.tree, format = 1)

print("Parsing dates")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))


for node in tree.traverse("postorder"):
	if node.is_leaf():
		if node.name[:4] == "node":
			node.detach()
		elif excess_leaf(node, from_date = options.from_date, to_date = options.to_date):
			node.detach()

tree.write(format=1, outfile=options.output)
