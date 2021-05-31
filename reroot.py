
import optparse
from ete3 import Tree
import pandas as pd

parser = optparse.OptionParser()
parser.add_option("--tree", dest='tree', type=str, help="full path to the newick file")
parser.add_option("--dates", dest='dates', type=str, help="full path to the csv file")
parser.add_option("--root", dest='root', type=str, help="name of the strain to use as an outgroup", default="")
parser.add_option("--output", dest='output', type=str, help="full path to the output newick file")
options, args = parser.parse_args()

t = Tree(options.tree, format=1) #flexible with internal node names
if options.root:
	if options.dates:
		raise Exception("You should either give me the root strain ( --root ), or a path to the file with dates ( --dates ), so i could find the oldest one; but not both!")
	t.set_outgroup( t&options.root )
else:
	if not options.dates:
		raise Exception("I need the root strain ( --root ) or the file with dates ( --dates ); recieved neither!") 
	else:
		print("Parsing dates")
		dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
		dates_dict = dict(zip(dates['seq_id'], dates['date']))
		tree_dates = {}
		for leaf in t.get_leaves():
			date = dates_dict.get(leaf.name, "unknown")
			if len(date)<10:
				date = "unknown"
			tree_dates[leaf.name] = date
		sorted_dates = sorted(tree_dates.items(), key=lambda item: item[1])
		print(sorted_dates[0:5])
		print(sorted_dates[0][0])
		t.set_outgroup( t&sorted_dates[0][0] )

	
	t.write(format=1, outfile=options.output)