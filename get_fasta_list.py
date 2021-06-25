import optparse

parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser
parser.add_option('-d', '--duplicates', help='file with duplicates (and cluster memebers)', type='str')
parser.add_option('-c', '--countries', help='file with ids and countries, output from meta_to_states.py', type='str')

options, args = parser.parse_args()

print("Parsing tree..")
tree = Tree(options.tree, format = 1)