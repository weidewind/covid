from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse


parser = argparse.ArgumentParser(description='Collect the results of global analysis into xlsx')
parser.add_argument("--fasta", dest='fasta', type=str, help="full path to the alignment file")
parser.add_argument("--output", dest='output', type=str, help="full path to the output file")
parser.add_argument("--ref_ind", dest='ref_ind', type=int, default=0, help="reference sequence index, default 0")
parser.add_argument("--start", dest='start', type=int, help="protein start")
parser.add_argument("--end", dest='end', type=int, help="protein end (position of the last letter)")
args = parser.parse_args()

# alignment ind -> site number in seq
def seq_numbering(seq):
	numbering = []
	num = 0
	for letter in list(refseq):
		if (letter != "-"):
			num += 1
		numbering.append(num)
	return(numbering)

align = AlignIO.read(args.fasta, "fasta")

refseq = align[args.ref_ind].seq
ref_numbering = seq_numbering(refseq)
prot = align[:, ref_numbering.index(args.start):ref_numbering.index(args.end+1)]

#AlignIO.write(prot, args.output, "fasta")
with open(args.output, "w") as out:
	for rec in prot:
		add = 45 - len(rec.id)
		out.write(rec.id + add*" ")
		out.write(str(rec.seq) + "\n")


