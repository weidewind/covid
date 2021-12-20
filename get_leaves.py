from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import optparse



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-o', '--output', help='', type='str')

options, args = parser.parse_args()
tree = Tree(options.tree, format=1)
with open(options.output, "w") as out:
	for l in tree.get_leaves():
		out.write(l.name + "\n")
