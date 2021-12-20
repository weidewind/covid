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
parser.add_option('-d', '--dates', help='file with ids and dates, output from meta_to_dates.py', type='str')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()
print(options)

selected_strains = []

meta_dict = get_country_dict()
print(dict(list(meta_dict.items())[0:5]))
pango = pd.read_csv(options.pangolined, sep=",")
print(pango[0:3])
pangolin_dict = dict(zip(pango['taxon'], pango['lineage']))
print(dict(list(pangolin_dict.items())[0:5]))

melted = pd.read_csv(options.melted_entries, sep=";")
print(melted[0:3])
melted_dict = dict(zip(melted['strain'], melted['entry']))
print(dict(list(melted_dict.items())[0:5]))

print("Parsing dates")
dates = pd.read_csv(options.dates, sep="\t", names=['seq_id', 'date'])
dates_dict = dict(zip(dates['seq_id'], dates['date']))
print(dict(list(dates_dict.items())[0:5]))

print("Parsing tree..")
tree = Tree(options.tree, format=1)
leaves = tree.get_leaves()


selected_set = []
for fasta in options.fasta.split(","):
	with open(fasta) as handle:
		print("Searching in " + fasta + "..")
		for record in SeqIO.parse(handle, "fasta"):
			if record.id in [l.name for l in leaves]:
				if meta_dict.get(record.id, "unknown") == "Russia" and (pangolin_dict.get(record.id, "unknown") == "B.1.617.2" or pangolin_dict.get(record.id, "unknown")[0:2] == "AY"):
					#out.write(record.id + " " + str(record.seq[1047:1048]) + " " + str(record.seq[27526:27527])+ "\n")
					#if record.seq[1047:1048].upper() == "T" and record.seq[27526:27527].upper() == "T":
					#if not record.seq[265:266].upper() == "A" or not record.seq[27393:27394].upper() == "A":
					if record.id not in selected_set:
						selected_set.append(record.id)
						selected_strains.append(record)



with open(options.output, "w") as out:
	out.write("taxon,n_1048,n_27527,is_double_mut,date,lineage,entry\n")
	for rec in selected_strains:
		n_1048 = str(rec.seq[1047:1048].upper())
		n_27527 = str(rec.seq[27526:27527].upper())
		double_mut = 1 if n_1048 == "T" and n_27527 == "T" else 0
		out.write(",".join([rec.id,n_1048,n_27527,str(double_mut),str(dates_dict.get(rec.id,"unknown")),pangolin_dict.get(rec.id, "unknown"),melted_dict.get(rec.id, "unknown")]) + "\n")

#SeqIO.write(mito_frags, "mitofrags.fasta", "fasta")

