#from Bio import Phylo
from ete3 import Tree
import time
import pandas as pd
import optparse
import sys
sys.setrecursionlimit(3500)


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-s','--steps', help='how many steps up', type='int', default=1)
parser.add_option('-e', '--entries_file', help='transmission_lineages.withduplicates.out.entries', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()



tree = Tree(options.tree, format=1)


with open(options.output, "w") as out:
    out.write("entry,node,distance\n")
    with open(options.entries_file, "r") as trf:
        for line in trf:
            entry_node = tree.search_nodes(name=line.strip())[0]
            dist = 0
            node = entry_node
            for i in range(0, options.steps):
                dist = dist + node.dist
                if not node.is_root():
                    node = node.up
                    out.write(",".join([entry_node.name, node.name, str(0-dist)]) + "\n")
