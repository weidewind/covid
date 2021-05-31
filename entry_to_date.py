import dendropy
from dendropy import TaxonSet, Tree, TreeList
from pprint import pprint
import optparse


parser = optparse.OptionParser()
parser.add_option('-d', '--tree_with_dates', help='nexus tree file with dates', type='str')
parser.add_option('-t', '--tree', help='newick tree file, topologically identical, but with other labels for internal nodes ', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


tree_dated = dendropy.Tree.get_from_path(options.tree_with_dates, "nexus", extract_comment_metadata=True, rooting='force-rooted')
tree_dated.write(path=options.output + ".tre", schema="newick")
tree = dendropy.Tree.get_from_path(options.tree, "newick", rooting='force-rooted')

with open(options.output, "w") as out:
	tree_iter = tree.postorder_node_iter()
	for nd in tree_dated.postorder_node_iter():
		#pprint(vars(n))
		n = next(tree_iter)
		if nd.is_leaf():
			out.write(n.taxon.label.replace(" ", "_") + "," + nd.taxon.label.replace(" ", "_") + "," + nd.annotations[0].value + "\n")
		else:
			out.write(n.label.replace(" ", "_") + "," + nd.label.replace(" ", "_") + "," + nd.annotations[0].value.replace(" ", "_") + "\n")
