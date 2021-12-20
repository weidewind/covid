import optparse
from ete3 import Tree
from meta import get_country_dict
import pandas as pd
from Bio import SeqIO


parser = optparse.OptionParser()
parser.add_option('-f', '--fasta', help='alignment file', type='str')
parser.add_option('-l', '--length', help='', type='int')
parser.add_option('-o', '--output', help='output', type='str')

options, args = parser.parse_args()

with open(options.fasta) as fasta:
	with open(options.output, "w") as out:
		for record in SeqIO.parse(fasta, "fasta"):
			out.write(">" + record.id + "\n")
			out.write(str(record.seq[0:options.length]) + "\n")
