from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import optparse



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-n', '--node', help='node', type='str')

options, args = parser.parse_args()
tree = Tree(options.tree, format=1)
print(str(tree.search_nodes(name=options.node)[0].up.name))
