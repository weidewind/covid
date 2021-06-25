from meta import get_date_dict
import datetime
from ete3 import Tree
import pandas as pd
import optparse
import re


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-m', '--meta', help='csv meta file form gisaid', type='str')
parser.add_option('-r', '--rpn_meta', help='xlsx rpn meta file', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format=1)

meta_dict = get_date_dict()

print("Writing output..")
with open(options.output, "w") as out:
	for strain, date in meta_dict.items():
		date = re.sub('-XX', '', date)
		if (len(date) < 10):
			date = "unknown"
		else:
			try:
				datetime.datetime.strptime(date, '%Y-%m-%d')
			except ValueError:
				try:
					datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S')
				except ValueError:
						print("Incorrect data format, should be YYYY-MM-DD [HH:MM:SS]: " + date)
						date = "unknown"
		out.write(strain + "\t" + date + "\n")
	for leaf in tree.iter_leaves():
		if leaf.name not in meta_dict:
			out.write(leaf.name + "\t" + "unknown" + "\n")
