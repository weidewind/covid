import optparse
from ete3 import Tree
from meta import get_country_dict
import pandas as pd
from Bio import SeqIO


parser = optparse.OptionParser()
parser.add_option('-f', '--fasta', help='msa', type='str')
parser.add_option('-d','--duplicates', type='str')
parser.add_option('-t', '--tree', type=str)
parser.add_option('-s', '--sites', help='sites, separated by comma', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()
print(options)

print("Parsing tree..")
tree = Tree(options.tree, format=1)
leaves = tree.get_leaves()

print("parsing duplicates..")
centroid_dict = {}
with open(options.duplicates) as dups:
	for line in dups:
		(centroid,duplicates) = line.strip().split("\t")
		for d in duplicates.split(";"):
			centroid_dict[d] = centroid


seq_dict = {}

for record in SeqIO.parse(options.fasta, "fasta"):
	seq = ""
	for s in options.sites.split(","):
		subs = str(record.seq[int(s)-1:int(s)]).upper()
		if subs not in ("A", "C", "G", "T"):
			subs  = "?"
		seq = seq + subs
	seq_dict[record.id] = seq

with open(options.output, "w") as out:
	for l in leaves:
		out.write(">" + l.name + "\n")
		if l.name in seq_dict:
			out.write(seq_dict[l.name] + "\n")
		else:
			out.write(seq_dict[centroid_dict[l.name]] + "\n")
