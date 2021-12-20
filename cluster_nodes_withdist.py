#from Bio import Phylo
from ete3 import Tree
import time
import pandas as pd
import optparse
import sys
sys.setrecursionlimit(3500)


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='tree file', type='str')
parser.add_option('-e', '--entries_file', help='transmission_lineages.withduplicates.out.entries', type='str')
parser.add_option('-x', '--exports_file', help='translin.entries.exports', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


def collect_distances_no_exports(entry_node, node, dist, dframe):
    print("entry " + entry_node.name + " node " + node.name)
    if entry_node.name not in e_to_exs or node.name not in e_to_exs[entry_node.name]:
        if not entry_node.name == node.name:
            dist = dist + node.dist
            print("init:")
            print(dframe[:min(dframe.shape[0],3)])
            dframe = dframe.append({"entry": entry_node.name, "node": node.name, "distance": dist}, ignore_index=True)
            print(dframe[:min(dframe.shape[0],3)])
            print("--------")
        if not node.is_leaf():
            for ch in node.get_children():
                dframe = collect_distances(entry_node, ch, dist, dframe)
            print("after processing children of  " + node.name + " returning ")
            return(dframe)
        else:
            print("after leaf " + node.name + " returning ")
            print(dframe[:min(dframe.shape[0],3)])
            return(dframe)
    else:
        print("after export node " + node.name + " returning ")
        print(dframe[:min(dframe.shape[0],3)])
        return(dframe)


def collect_distances(entry_node, node, dist, dframe):
    print("entry " + entry_node.name + " node " + node.name)
    if not entry_node.name == node.name:
        ntype = "internal" if (entry_node.name not in e_to_exs or node.name not in e_to_exs[entry_node.name]) else "export" 
        dist = dist + node.dist
        print("init:")
        print(dframe[:min(dframe.shape[0],3)])
        dframe = dframe.append({"entry": entry_node.name, "node": node.name, "distance": dist, "type": ntype}, ignore_index=True)
        print(dframe[:min(dframe.shape[0],3)])
        print("--------")
    if not node.is_leaf() and (entry_node.name not in e_to_exs or node.name not in e_to_exs[entry_node.name]):
        for ch in node.get_children():
            dframe = collect_distances(entry_node, ch, dist, dframe)
        print("after processing children of  " + node.name + " returning ")
        return(dframe)
    else:
        print("after leaf or export " + node.name + " returning ")
        print(dframe[:min(dframe.shape[0],3)])
        return(dframe)


tree = Tree(options.tree, format=1)

e_to_exs = {}
with open(options.exports_file, "r") as exf:
    for line in exf:
        (entry, ex) = line.strip().split("\t")
        if entry not in e_to_exs:
            e_to_exs[entry] = []
        e_to_exs[entry].append(ex)

df = pd.DataFrame(columns=["entry", "node", "distance"])
with open(options.entries_file, "r") as trf:
    for line in trf:
        entry_node = tree.search_nodes(name=line.strip())[0]
        df = collect_distances(entry_node, entry_node, dist=0, dframe=df)
df.to_csv(options.output, index=False)

