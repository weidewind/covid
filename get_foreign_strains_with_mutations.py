import optparse
from ete3 import Tree
from meta import get_country_dict
import pandas as pd
from Bio import SeqIO


parser = optparse.OptionParser()
parser.add_option('-t', '--tree', help='newick', type='str')
parser.add_option('-f', '--fasta', help='alignment files, separated by comma', type='str')
parser.add_option('-p', '--pangolined', help='', type='str')
parser.add_option('-m', '--melted_entries', help='', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()


pango = pd.read_csv(options.pangolined, sep=",")
print(pango[0:3])
pangolin_dict = dict(zip(pango['taxon'], pango['lineage']))
print(dict(list(pangolin_dict.items())[0:5]))

melted = pd.read_csv(options.melted_entries, sep=";")
print(melted[0:3])

print("Parsing tree..")
tree = Tree(options.tree, format=1)
leaves = tree.get_leaves()

strain_dict = {}
for fasta in options.fasta.split(","):
	with open(fasta) as handle:
		print("Searching in " + fasta + "..")
		for record in SeqIO.parse(handle, "fasta"):
			if record.id in [l.name for l in leaves]:
				strain_dict[record.id] = record



with open(options.output, "w") as out:
	out.write("taxon,n_1048,n_27527,is_double_mut,date,lineage,entry\n")
	for index, row in melted.iterrows():
		rec = strain_dict[row["strain"]]
		n_1048 = str(rec.seq[1047:1048].upper())
		n_27527 = str(rec.seq[27526:27527].upper())
		double_mut = 1 if n_1048 == "T" and n_27527 == "T" else 0
		out.write(",".join([row["strain"],n_1048,n_27527,str(double_mut),row["date"],pangolin_dict.get(rec.id, "unknown"),row["entry"]]) + "\n")

#SeqIO.write(mito_frags, "mitofrags.fasta", "fasta")

