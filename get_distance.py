from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import optparse



parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-n', '--nodes', help='nodes separated by comma', type='str')

options, args = parser.parse_args()
tree = Tree(options.tree, format=1)
(n1,n2) = options.nodes.split(",")
print(str(tree.get_distance(n1,n2)))
